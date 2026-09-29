#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


MODES = [
    ("M000", False, False, False),
    ("M001", False, False, True),
    ("M010", False, True,  False),
    ("M011", False, True,  True),
    ("M100", True,  False, False),
    ("M101", True,  False, True),
    ("M110", True,  True,  False),
    ("M111", True,  True,  True),
]


def as_bool(x):
    return str(x).strip() in {"1", "true", "True", "TRUE", "yes", "Yes"}


def role_key(primary, secondary, supplementary):
    return (bool(primary), bool(secondary), bool(supplementary))


def record_key(qname, start, primary, secondary, supplementary):
    return (qname, int(start), *role_key(primary, secondary, supplementary))


def segment_tuple_from_expected(row):
    return (
        row["expected_type"],
        row["expected_op"],
        int(row["expected_length"]),
        int(row["expected_ref_consumed"]),
        int(row["expected_read_consumed"]),
    )


def segment_tuple_from_observed(seg):
    return (
        seg.type,
        seg.op,
        int(seg.length),
        int(seg.ref_consumed),
        int(seg.read_consumed),
    )


def format_segments(segs):
    return ";".join(
        f"{typ}:{op}:{length}:ref={refc}:read={readc}"
        for typ, op, length, refc, readc in segs
    )


def expected_source(use_ref, use_md, use_cs, has_md, has_cs):
    # This selects which manually predefined truth column/source to compare.
    # It does NOT generate expected segments.
    if use_ref:
        return "ref"
    if use_cs and has_cs:
        return "cs"
    if use_md and has_md:
        return "md"
    return "cigar"


def compare_metadata(obs, exp):
    failures = []

    checks = [
        ("chrom", obs.chrom, exp["chrom"]),
        ("start", int(obs.start), int(exp["expected_start_0based"])),
        ("end", int(obs.end), int(exp["expected_end_0based"])),
        ("strand", "-" if obs.reverse else "+", exp["expected_strand"]),
        ("primary", bool(obs.primary), as_bool(exp["expected_primary"])),
        ("secondary", bool(obs.secondary), as_bool(exp["expected_secondary"])),
        ("supplementary", bool(obs.supplementary), as_bool(exp["expected_supplementary"])),
    ]

    for field, observed, expected in checks:
        if observed != expected:
            failures.append(f"{field}: expected={expected!r}, observed={observed!r}")

    return failures


def compare_segments(observed_segments, expected_segments):
    failures = []

    if len(observed_segments) != len(expected_segments):
        failures.append(
            f"segment_count: expected={len(expected_segments)}, observed={len(observed_segments)}"
        )

    for i, (exp, obs) in enumerate(zip(expected_segments, observed_segments), start=1):
        if exp != obs:
            failures.append(
                f"segment_{i}: expected={exp!r}, observed={obs!r}"
            )

    return failures


def main():
    parser = argparse.ArgumentParser(
        description="Validate BamSnap-LRS parsing against manually predefined synthetic truth."
    )
    parser.add_argument("--bam", required=True, help="Synthetic parsing BAM")
    parser.add_argument("--fa", required=True, help="Synthetic reference FASTA")
    parser.add_argument("--records-truth", required=True, help="parsing_expected_records.tsv")
    parser.add_argument("--segments-truth", required=True, help="parsing_expected_segments.tsv")
    parser.add_argument(
        "--bamsnap-src",
        default=None,
        help="Path to Bamsnap-LRS/src. If omitted, bamsnap_lrs must already be importable."
    )
    parser.add_argument("--chrom", default="chrTest", help="Chromosome to fetch (default: chrTest)")
    parser.add_argument("--start", type=int, default=0, help="Fetch start, 0-based (default: 0)")
    parser.add_argument("--end", type=int, default=1000, help="Fetch end, 0-based half-open (default: 1000)")
    parser.add_argument(
        "-o", "--output",
        default="parsing_validation_results.tsv",
        help="Output comparison TSV"
    )
    args = parser.parse_args()

    if args.bamsnap_src:
        sys.path.insert(0, str(Path(args.bamsnap_src).resolve()))

    try:
        from bamsnap_lrs.reader import fetch_reads
    except ImportError as e:
        raise SystemExit(
            "Cannot import bamsnap_lrs.reader. Use --bamsnap-src /path/to/Bamsnap-LRS/src "
            "or set PYTHONPATH before running.\n"
            f"Original error: {e}"
        )

    # ----------------------------
    # Load record-level manual truth
    # ----------------------------
    with open(args.records_truth, newline="") as f:
        record_rows = list(csv.DictReader(f, delimiter="\t"))

    records_by_id = {}
    truth_key_to_id = {}

    for row in record_rows:
        rid = row["record_id"]
        if rid in records_by_id:
            raise ValueError(f"Duplicate record_id in record truth: {rid}")

        key = record_key(
            row["qname"],
            row["expected_start_0based"],
            as_bool(row["expected_primary"]),
            as_bool(row["expected_secondary"]),
            as_bool(row["expected_supplementary"]),
        )
        if key in truth_key_to_id:
            raise ValueError(
                f"Record truth is not uniquely identifiable: {rid} and {truth_key_to_id[key]}"
            )

        records_by_id[rid] = row
        truth_key_to_id[key] = rid

    # ----------------------------
    # Load segment-level manual truth
    # ----------------------------
    expected_segments = defaultdict(list)

    with open(args.segments_truth, newline="") as f:
        segment_rows = list(csv.DictReader(f, delimiter="\t"))

    for row in segment_rows:
        key = (row["record_id"], row["source"])
        expected_segments[key].append(
            (int(row["segment_index"]), segment_tuple_from_expected(row))
        )

    # Sort by manually assigned segment_index.
    expected_segments = {
        key: [seg for _, seg in sorted(items, key=lambda x: x[0])]
        for key, items in expected_segments.items()
    }

    results = []

    # ----------------------------
    # Run all 8 parsing configurations
    # ----------------------------
    for mode, use_ref, use_md, use_cs in MODES:
        reads = fetch_reads(
            bam_path=args.bam,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
            mapq_min=0,
            show_supp=True,
            show_secondary=True,
            use_md=use_md,
            use_cs=use_cs,
            use_ref=use_ref,
            fa_path=args.fa,
        )

        observed_by_id = {}
        unmatched_observed = []

        for read in reads:
            key = record_key(
                read.qname,
                read.start,
                read.primary,
                read.secondary,
                read.supplementary,
            )
            rid = truth_key_to_id.get(key)
            if rid is None:
                unmatched_observed.append(
                    f"{read.qname}@{read.start}"
                    f"(primary={read.primary},secondary={read.secondary},supp={read.supplementary})"
                )
                continue
            if rid in observed_by_id:
                raise ValueError(f"{mode}: duplicate observed record matched to {rid}")
            observed_by_id[rid] = read

        # One result row per expected BAM alignment record.
        for rid, exp in records_by_id.items():
            failures = []
            read = observed_by_id.get(rid)

            has_md = as_bool(exp["has_md"])
            has_cs = as_bool(exp["has_cs"])
            source = expected_source(use_ref, use_md, use_cs, has_md, has_cs)

            exp_segs = expected_segments.get((rid, source))
            if exp_segs is None:
                failures.append(f"missing_manual_truth_for_source:{source}")
                exp_segs = []

            if read is None:
                failures.append("record_missing_from_fetch_reads")
                obs_segs = []
                observed_start = ""
                observed_end = ""
                observed_strand = ""
                observed_primary = ""
                observed_secondary = ""
                observed_supplementary = ""
            else:
                failures.extend(compare_metadata(read, exp))
                obs_segs = [segment_tuple_from_observed(s) for s in read.segments]
                failures.extend(compare_segments(obs_segs, exp_segs))

                observed_start = read.start
                observed_end = read.end
                observed_strand = "-" if read.reverse else "+"
                observed_primary = int(bool(read.primary))
                observed_secondary = int(bool(read.secondary))
                observed_supplementary = int(bool(read.supplementary))

            results.append({
                "record_id": rid,
                "qname": exp["qname"],
                "mode": mode,
                "use_ref": int(use_ref),
                "use_md": int(use_md),
                "use_cs": int(use_cs),
                "has_md": int(has_md),
                "has_cs": int(has_cs),
                "expected_source": source,

                "expected_start": exp["expected_start_0based"],
                "observed_start": observed_start,
                "expected_end": exp["expected_end_0based"],
                "observed_end": observed_end,
                "expected_strand": exp["expected_strand"],
                "observed_strand": observed_strand,
                "expected_primary": exp["expected_primary"],
                "observed_primary": observed_primary,
                "expected_secondary": exp["expected_secondary"],
                "observed_secondary": observed_secondary,
                "expected_supplementary": exp["expected_supplementary"],
                "observed_supplementary": observed_supplementary,

                "expected_segments": format_segments(exp_segs),
                "observed_segments": format_segments(obs_segs),

                "status": "PASS" if not failures else "FAIL",
                "failure_reason": " | ".join(failures),
            })

        # Preserve unexpected records as explicit failures.
        for desc in unmatched_observed:
            results.append({
                "record_id": "",
                "qname": "",
                "mode": mode,
                "use_ref": int(use_ref),
                "use_md": int(use_md),
                "use_cs": int(use_cs),
                "has_md": "",
                "has_cs": "",
                "expected_source": "",
                "expected_start": "",
                "observed_start": "",
                "expected_end": "",
                "observed_end": "",
                "expected_strand": "",
                "observed_strand": "",
                "expected_primary": "",
                "observed_primary": "",
                "expected_secondary": "",
                "observed_secondary": "",
                "expected_supplementary": "",
                "observed_supplementary": "",
                "expected_segments": "",
                "observed_segments": "",
                "status": "FAIL",
                "failure_reason": f"unexpected_record:{desc}",
            })

    # ----------------------------
    # Write machine-readable results
    # ----------------------------
    fieldnames = [
        "record_id","qname","mode","use_ref","use_md","use_cs",
        "has_md","has_cs","expected_source",
        "expected_start","observed_start","expected_end","observed_end",
        "expected_strand","observed_strand",
        "expected_primary","observed_primary",
        "expected_secondary","observed_secondary",
        "expected_supplementary","observed_supplementary",
        "expected_segments","observed_segments",
        "status","failure_reason"
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    passed = sum(r["status"] == "PASS" for r in results)
    failed = sum(r["status"] == "FAIL" for r in results)

    print(f"Output: {args.output}")
    print(f"Expected records: {len(records_by_id)}")
    print(f"Parsing modes: {len(MODES)}")
    print(f"Expected record × mode tests: {len(records_by_id) * len(MODES)}")
    print(f"Result rows: {len(results)}")
    print(f"PASS: {passed}")
    print(f"FAIL: {failed}")

    if failed:
        print("\nFailed tests:")
        for r in results:
            if r["status"] == "FAIL":
                rid = r["record_id"] or "<unexpected>"
                print(f"  {rid}\t{r['mode']}\t{r['failure_reason']}")
        sys.exit(1)


if __name__ == "__main__":
    main()
