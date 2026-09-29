#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
import pysam

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MISMATCH_RE = re.compile(r"offset\s+(\d+):\s*([ACGT])\s*=>\s*([ACGT])", re.I)
REF_RANGE_RE = re.compile(r"ref\[(\d+):(\d+)\)")
INSERT_RE = re.compile(r"\binsert\s+([ACGT]+)\b", re.I)
SOFTCLIP_RE = re.compile(r"soft clips\s+([ACGTN]+)/([ACGTN]+)", re.I)

def read_fasta(path):
    references = {}
    name = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    references[name] = "".join(seq)

                name = line[1:].split()[0]
                if not name:
                    raise ValueError("Empty FASTA sequence name.")
                if name in references:
                    raise ValueError(f"Duplicate FASTA sequence name: {name}")
                seq = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before the first header.")
                seq.append(line.upper())

    if name is not None:
        references[name] = "".join(seq)

    if not references:
        raise ValueError("No FASTA sequence found.")

    return references

def parse_cigar(cigar):
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise ValueError(f"Unsupported or malformed CIGAR: {cigar}")
    return ops

def query_length_from_cigar(cigar):
    return sum(n for n, op in parse_cigar(cigar) if op in "MIS=X")

def ref_length_from_cigar(cigar):
    return sum(n for n, op in parse_cigar(cigar) if op in "MDN=X")

def parse_mismatches(event):
    out = {}
    for off, ref_base, query_base in MISMATCH_RE.findall(event):
        out[int(off)] = (ref_base.upper(), query_base.upper())
    return out

def parse_insertions(event):
    return [x.upper() for x in INSERT_RE.findall(event)]

def parse_softclips(event):
    m = SOFTCLIP_RE.search(event)
    return [m.group(1).upper(), m.group(2).upper()] if m else []

def build_reference_concat_query(event, ref):
    ranges = REF_RANGE_RE.findall(event)
    if not ranges:
        return None
    return "".join(ref[int(s):int(e)] for s, e in ranges)

def default_softclip(index, length):
    # R07 specifies clip lengths but not clip bases. These bases are not part
    # of the clipping truth being tested, so use deterministic filler.
    bases = ("A", "T", "C", "G")
    return bases[index % len(bases)] * length

def build_cigar_derived_query(row, ref):
    start = int(row["ref_start_0based"])
    cigar = row["cigar"]
    event = row["manual_event(0-based offset)"]

    mismatches = parse_mismatches(event)
    insertions = parse_insertions(event)
    softclips = parse_softclips(event)

    ref_pos = start
    aligned_offset = 0
    ins_i = 0
    soft_i = 0
    q = []

    for n, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            for _ in range(n):
                ref_base = ref[ref_pos]
                query_base = ref_base

                if aligned_offset in mismatches:
                    expected_ref, query_base = mismatches[aligned_offset]
                    if ref_base != expected_ref:
                        raise ValueError(
                            f'{row["record_id"]}: mismatch truth disagrees with reference at '
                            f'0-based position {ref_pos}: table says {expected_ref}, FASTA has {ref_base}'
                        )
                    if query_base == ref_base:
                        raise ValueError(f'{row["record_id"]}: mismatch query base equals reference base.')

                if op == "X" and aligned_offset not in mismatches:
                    raise ValueError(
                        f'{row["record_id"]}: CIGAR X at aligned offset {aligned_offset} '
                        f'has no explicit REF=>QUERY event in the table.'
                    )
                if op == "=" and aligned_offset in mismatches:
                    raise ValueError(
                        f'{row["record_id"]}: mismatch event falls inside a CIGAR "=" block.'
                    )

                q.append(query_base)
                ref_pos += 1
                aligned_offset += 1

        elif op == "I":
            if ins_i >= len(insertions):
                raise ValueError(f'{row["record_id"]}: CIGAR contains {n}I but no insertion sequence is defined.')
            ins = insertions[ins_i]
            if len(ins) != n:
                raise ValueError(
                    f'{row["record_id"]}: insertion sequence {ins} has length {len(ins)}, expected {n}.'
                )
            q.append(ins)
            ins_i += 1

        elif op in ("D", "N"):
            ref_pos += n
            aligned_offset += n

        elif op == "S":
            if soft_i < len(softclips):
                clip = softclips[soft_i]
                if len(clip) != n:
                    raise ValueError(
                        f'{row["record_id"]}: soft clip sequence {clip} has length {len(clip)}, expected {n}.'
                    )
            else:
                clip = default_softclip(soft_i, n)
            q.append(clip)
            soft_i += 1

        elif op in ("H", "P"):
            pass

    query = "".join(q)
    expected_len = query_length_from_cigar(cigar)
    if len(query) != expected_len:
        raise ValueError(
            f'{row["record_id"]}: generated query length {len(query)} != CIGAR query length {expected_len}.'
        )
    return query

def validate_manual_reference_events(row, ref):
    event = row["manual_event(0-based offset)"]

    # Validate explicitly written ref[start:end)=SEQ statements.
    for s, e in REF_RANGE_RE.findall(event):
        token = f"ref[{s}:{e})="
        if token in event:
            after = event.split(token, 1)[1]
            expected = re.match(r"([ACGT]+)", after, re.I)
            if expected:
                observed = ref[int(s):int(e)]
                if observed != expected.group(1).upper():
                    raise ValueError(
                        f'{row["record_id"]}: table says ref[{s}:{e})={expected.group(1).upper()}, '
                        f'but FASTA has {observed}.'
                    )

def validate_tag_design(row):
    md = row["manual_MD"].strip()
    cs = row["manual_cs"].strip()
    design = row["tag_design"].strip().lower()

    expected = {
        "md+cs": (True, True),
        "md only": (True, False),
        "cs only": (False, True),
        "neither": (False, False),
    }
    if design not in expected:
        raise ValueError(f'{row["record_id"]}: unknown tag_design: {row["tag_design"]}')

    want_md, want_cs = expected[design]
    if bool(md) != want_md or bool(cs) != want_cs:
        raise ValueError(
            f'{row["record_id"]}: tag_design={row["tag_design"]} is inconsistent with manual_MD/manual_cs.'
        )

def main():
    parser = argparse.ArgumentParser(
        description="Generate the BamSnap-LRS synthetic parsing-validation BAM from a manually specified TSV."
    )
    parser.add_argument("-r", "--reference", required=True, help="Synthetic reference FASTA")
    parser.add_argument("-d", "--design", required=True, help="Manual CIGAR/MD/cs design TSV")
    parser.add_argument("-o", "--output", default="synthetic_parsing.bam", help="Output coordinate-sorted BAM")
    parser.add_argument("--sam", default=None, help="Optional SAM copy for manual inspection")
    args = parser.parse_args()

    references = read_fasta(args.reference)
    ref_names = list(references)
    ref_id_by_name = {name: i for i, name in enumerate(ref_names)}

    with open(args.design, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    required = {
        "case_id", "record_id", "qname", "flag", "chrom_name", "ref_start_0based", "mapq", "cigar",
        "manual_event(0-based offset)", "manual_MD", "manual_cs", "tag_design", "note"
    }
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Missing TSV columns: {sorted(missing)}")

    # Build any explicitly defined shared molecule sequence first.
    shared_query = {}
    for row in rows:
        chrom = row["chrom_name"]
        if chrom not in references:
            raise ValueError(
                f'{row["record_id"]}: chromosome {chrom!r} is not present in the input FASTA.'
            )
        ref = references[chrom]
        event = row["manual_event(0-based offset)"]
        q = build_reference_concat_query(event, ref)
        if q is not None and event.strip().startswith("query="):
            shared_query[row["qname"]] = q

    records = []
    for row in rows:
        validate_tag_design(row)

        flag = int(row["flag"])
        chrom = row["chrom_name"]
        if chrom not in references:
            raise ValueError(
                f'{row["record_id"]}: chromosome {chrom!r} is not present in the input FASTA.'
            )
        ref = references[chrom]
        ref_id = ref_id_by_name[chrom]
        validate_manual_reference_events(row, ref)
        start = int(row["ref_start_0based"])
        mapq = int(row["mapq"])
        cigar = row["cigar"]
        event = row["manual_event(0-based offset)"]

        if "same query" in event.lower():
            if row["qname"] not in shared_query:
                raise ValueError(f'{row["record_id"]}: no shared query recipe found for QNAME {row["qname"]}.')
            query = shared_query[row["qname"]]
        elif row["qname"] in shared_query:
            query = shared_query[row["qname"]]
        else:
            query = build_cigar_derived_query(row, ref)

        if start < 0 or start + ref_length_from_cigar(cigar) > len(ref):
            raise ValueError(f'{row["record_id"]}: alignment extends outside the reference.')

        if len(query) != query_length_from_cigar(cigar):
            raise ValueError(
                f'{row["record_id"]}: query length {len(query)} does not match CIGAR query length '
                f'{query_length_from_cigar(cigar)}.'
            )

        a = pysam.AlignedSegment()
        a.query_name = row["qname"]
        a.flag = flag
        a.reference_id = ref_id
        a.reference_start = start
        a.mapping_quality = mapq
        a.cigarstring = cigar
        a.query_sequence = query
        a.query_qualities = pysam.qualitystring_to_array("I" * len(query))

        md = row["manual_MD"].strip()
        cs = row["manual_cs"].strip()
        if md:
            a.set_tag("MD", md, value_type="Z")
        if cs:
            a.set_tag("cs", cs, value_type="Z")

        records.append((ref_id, start, flag, row["record_id"], a))

    # Coordinate-sort before indexing.
    records.sort(key=lambda x: (x[0], x[1], x[2] & 0x100 != 0, x[2] & 0x800 != 0, x[3]))

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(references[name])} for name in ref_names],
    }

    output = Path(args.output)
    with pysam.AlignmentFile(str(output), "wb", header=header) as bam:
        for _, _, _, _, a in records:
            bam.write(a)

    pysam.index(str(output))

    if args.sam:
        with pysam.AlignmentFile(str(output), "rb") as bam, pysam.AlignmentFile(
            args.sam, "w", header=bam.header
        ) as sam:
            for rec in bam:
                sam.write(rec)

    print(f"Generated: {output}")
    print(f"Indexed:   {output}.bai")
    if args.sam:
        print(f"SAM copy:  {args.sam}")
    print(f"Records:   {len(records)}")

if __name__ == "__main__":
    main()
