#!/usr/bin/env python3
"""
Validate one BamSnap-LRS intermediate TSV against two precomputed reference files:

1) SAM-reference TSV made from `samtools view` output.
   Accepted formats (tab-delimited, no header required):
       QNAME  FLAG  RNAME  POS  MAPQ
   or, preferably:
       QNAME  FLAG  RNAME  POS  MAPQ  CIGAR
   POS is SAM 1-based and is converted internally to 0-based start.

2) PAF made from the same regional SAM records with paftools.js sam2paf.
   cg:Z: is used when available to compare alignment-event semantics.

3) BamSnap-LRS compact intermediate TSV.

The script writes ONE summary TSV row to -o/--out.  It does not call samtools
and does not read BAM directly. Reference records absent from the BamSnap-LRS
intermediate output are reported as excluded records. Downstream grouping, order,
and visual-semantic checks are performed on the retained/matched record subset.
"""

import argparse
import csv
import sys
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
csv.field_size_limit(sys.maxsize)

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class SamRefRec:
    qname: str
    flag: int
    chrom: str
    start: int       # 0-based
    end: Optional[int]
    mapq: int
    cigar: Optional[str]
    strand: str
    primary: bool
    supplementary: bool
    secondary: bool


@dataclass
class PafRec:
    qname: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    chrom: str
    start: int
    end: int
    mapq: int
    cg: Optional[str]


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_optional_float(value) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip()
    if s in {"", ".", "None", "null", "NA", "nan"}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def cigar_ref_consumed(cigar: str) -> int:
    if not cigar or cigar == "*":
        return 0
    return sum(
        int(n)
        for n, op in _CIGAR_RE.findall(cigar)
        if op in {"M", "D", "N", "=", "X"}
    )


def read_sam_reference(path: str) -> Tuple[List[SamRefRec], bool]:
    """Read 5- or 6-column SAM-derived reference TSV.

    Returns (records, has_cigar_for_all_records).
    """
    out: List[SamRefRec] = []
    all_have_cigar = True
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                raise SystemExit(
                    f"Invalid SAM-reference line {line_no}: expected at least 5 tab-delimited fields"
                )
            # Optional header support.
            if line_no == 1 and f[0].upper() in {"QNAME", "READ", "READ_NAME"}:
                continue
            try:
                flag = int(f[1])
                pos1 = int(f[3])
                mapq = int(f[4])
            except ValueError as exc:
                raise SystemExit(
                    f"Invalid FLAG/POS/MAPQ on SAM-reference line {line_no}: {line.rstrip()}"
                ) from exc

            chrom = f[2]
            cigar = f[5] if len(f) >= 6 and f[5] not in {"", "."} else None
            if cigar is None:
                all_have_cigar = False
            start = pos1 - 1
            end = start + cigar_ref_consumed(cigar) if cigar else None
            supp = bool(flag & 0x800)
            sec = bool(flag & 0x100)
            prim = not supp and not sec
            out.append(SamRefRec(
                qname=f[0],
                flag=flag,
                chrom=chrom,
                start=start,
                end=end,
                mapq=mapq,
                cigar=cigar,
                strand="-" if (flag & 0x10) else "+",
                primary=prim,
                supplementary=supp,
                secondary=sec,
            ))
    return out, all_have_cigar


def parse_paf_tags(fields: Sequence[str]) -> Dict[str, str]:
    tags = {}
    for x in fields:
        p = x.split(":", 2)
        if len(p) == 3:
            tags[p[0]] = p[2]
    return tags


def read_paf(path: str) -> List[PafRec]:
    out: List[PafRec] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise SystemExit(f"Invalid PAF line {line_no}: fewer than 12 fields")
            tags = parse_paf_tags(f[12:])
            out.append(PafRec(
                qname=f[0],
                query_length=int(f[1]),
                query_start=int(f[2]),
                query_end=int(f[3]),
                strand=f[4],
                chrom=f[5],
                start=int(f[7]),
                end=int(f[8]),
                mapq=int(f[11]),
                cg=tags.get("cg"),
            ))
    return out


def read_intermediate(path: str) -> List[dict]:
    required = {
        "qname", "chrom", "start", "end", "strand", "mapq",
        "primary", "supplementary", "secondary",
        "query_start", "query_end", "query_length",
        "layout_query_start", "layout_query_end",
        "qname_group_size", "query_order_in_group", "stack_row",
        "supplementary_group_color", "supplementary_group_opacity",
        "read_body_color", "read_body_opacity", "arrow_direction",
        "render_semantics_json",
    }
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Intermediate TSV is missing required columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            try:
                row["_render_semantics"] = json.loads(row["render_semantics_json"] or "[]")
            except json.JSONDecodeError as exc:
                raise SystemExit(
                    f"Invalid render_semantics_json for {row.get('qname', '?')}: {exc}"
                ) from exc
            rows.append(row)
    return rows


# -----------------------------------------------------------------------------
# Matching keys
# -----------------------------------------------------------------------------

def ref_key(rec: SamRefRec, strong: bool):
    if strong and rec.end is not None:
        return (rec.qname, rec.chrom, rec.start, rec.end, rec.strand, rec.mapq)
    return (rec.qname, rec.chrom, rec.start, rec.strand, rec.mapq)


def paf_key(rec: PafRec):
    return (rec.qname, rec.chrom, rec.start, rec.end, rec.strand, rec.mapq)


def tsv_strong_key(row: dict):
    return (
        row["qname"], row["chrom"], int(row["start"]), int(row["end"]),
        row["strand"], int(row["mapq"]),
    )


def tsv_weak_key(row: dict):
    return (
        row["qname"], row["chrom"], int(row["start"]),
        row["strand"], int(row["mapq"]),
    )


def multiset_match_count(left, right) -> int:
    return sum((Counter(left) & Counter(right)).values())


def occurrence_map(records, key_func):
    seen = Counter()
    out = {}
    for rec in records:
        k = key_func(rec)
        occ = seen[k]
        seen[k] += 1
        out[k + (occ,)] = rec
    return out


# -----------------------------------------------------------------------------
# Render semantics
# -----------------------------------------------------------------------------

def merge_events(events: List[Tuple[str, int]]) -> List[Tuple[str, int]]:
    out: List[Tuple[str, int]] = []
    for t, n in events:
        n = int(n)
        if out and out[-1][0] == t:
            out[-1] = (t, out[-1][1] + n)
        else:
            out.append((t, n))
    return out


def paf_cg_events(cg: Optional[str]):
    if not cg:
        return [], False
    mapping = {
        "=": "match", "X": "mismatch", "M": "generic_match",
        "I": "ins", "D": "del", "N": "ref_skip",
        "S": "soft", "H": "hard",
    }
    events = []
    detailed = False
    has_generic_m = False
    for n, op in _CIGAR_RE.findall(cg):
        if op in {"=", "X"}:
            detailed = True
        if op == "M":
            has_generic_m = True
        t = mapping.get(op)
        if t is None or t in {"soft", "hard"}:
            continue
        events.append((t, int(n)))
    return merge_events(events), detailed and not has_generic_m


def semantic_events(row: dict):
    items = [
        x for x in row["_render_semantics"]
        if x.get("type") != "supplementary_overview_body"
    ]
    events = []
    for x in items:
        if x.get("type") is None or x.get("length") is None:
            continue
        events.append((str(x["type"]), int(x["length"])))
    return merge_events(events), items


def genericize(events: List[Tuple[str, int]]) -> List[Tuple[str, int]]:
    out = []
    for t, n in events:
        if t in {"match", "mismatch", "generic_match"}:
            t = "aligned"
        out.append((t, n))
    return merge_events(out)


def compare_render_events(paf: PafRec, row: dict):
    expected, detailed = paf_cg_events(paf.cg)
    observed, items = semantic_events(row)
    if not expected:
        return None, 0, 0, 0, 0, 0, 0

    exp_cmp = expected if detailed else genericize(expected)
    obs_cmp = observed if detailed else genericize(observed)
    sequence_match = exp_cmp == obs_cmp
    event_tested = max(len(exp_cmp), len(obs_cmp))
    event_matching = sum(1 for a, b in zip(exp_cmp, obs_cmp) if a == b)

    symbol_tested = 0
    symbol_matching = 0
    label_tested = 0
    label_matching = 0
    for item in items:
        if not parse_bool(item.get("drawn", True)):
            continue
        t = item.get("type")
        if t in {"match", "mismatch", "ins", "del", "ref_skip"}:
            expected_symbol = "line" if t == "ref_skip" else "rect"
            symbol_tested += 1
            symbol_matching += int(item.get("symbol") == expected_symbol)
        if t in {"ins", "del"} and item.get("length") is not None:
            n = int(item["length"])
            expected_label = f"I({n})" if t == "ins" else str(n)
            label_tested += 1
            label_matching += int(item.get("label") == expected_label)

    return (
        sequence_match,
        event_tested, event_matching,
        symbol_tested, symbol_matching,
        label_tested, label_matching,
    )


# -----------------------------------------------------------------------------
# Group/layout/visual checks
# -----------------------------------------------------------------------------

def intervals_overlap(a_start, a_end, b_start, b_end) -> bool:
    return a_start < b_end and b_start < a_end


def expected_opacities(n: int) -> List[float]:
    if n <= 1:
        return [1.0]
    if n == 2:
        return [1.0, 0.6]
    return [1.0 - 0.8 * i / (n - 1) for i in range(n)]


def float_equal(a, b, tol=1e-6) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def main():
    ap = argparse.ArgumentParser(
        description="Validate BamSnap-LRS intermediate TSV against precomputed SAM-reference TSV and PAF."
    )
    ap.add_argument("--reference", required=True,
                    help="SAM-derived TSV: QNAME FLAG RNAME POS MAPQ [CIGAR]")
    ap.add_argument("--paf", required=True)
    ap.add_argument("--intermediate", required=True)
    ap.add_argument("-o", "--out", required=True, help="Output one-row validation summary TSV")
    args = ap.parse_args()

    refs, strong_reference = read_sam_reference(args.reference)
    pafs = read_paf(args.paf)
    rows = read_intermediate(args.intermediate)

    # Selection / identity counts.
    ref_keys = [ref_key(r, strong_reference) for r in refs]
    row_ref_keys = [tsv_strong_key(r) if strong_reference else tsv_weak_key(r) for r in rows]
    matched_ref = multiset_match_count(ref_keys, row_ref_keys)
    excluded_ref = len(refs) - matched_ref
    unexpected_ref = len(rows) - matched_ref

    paf_keys = [paf_key(p) for p in pafs]
    row_paf_keys = [tsv_strong_key(r) for r in rows]
    matched_paf = multiset_match_count(paf_keys, row_paf_keys)

    # Occurrence-aware reference-to-intermediate matching.
    ref_map = occurrence_map(refs, lambda r: ref_key(r, strong_reference))
    row_ref_map = occurrence_map(rows, tsv_strong_key if strong_reference else tsv_weak_key)

    flag_tested = flag_matches = 0
    arrow_tested = arrow_matches = 0
    end_tested = end_matches = 0
    failure_examples = []
    retained_refs = []
    retained_rows = []

    for k, ref in ref_map.items():
        row = row_ref_map.get(k)
        if row is None:
            # Keep the exclusion note concise. MAPQ=0 is reported explicitly
            # because it is a common pre-render exclusion category in these tests.
            if len(failure_examples) < 5:
                reason = "mapq0" if ref.mapq == 0 else "excluded"
                failure_examples.append(f"{reason}:{ref.qname}:{ref.chrom}:{ref.start}")
            continue

        retained_refs.append(ref)
        retained_rows.append(row)

        flag_tested += 1
        status_ok = (
            parse_bool(row["primary"]) == ref.primary
            and parse_bool(row["supplementary"]) == ref.supplementary
            and parse_bool(row["secondary"]) == ref.secondary
        )
        flag_matches += int(status_ok)
        if not status_ok and len(failure_examples) < 5:
            failure_examples.append(f"flag:{ref.qname}:{ref.chrom}:{ref.start}")

        arrow_tested += 1
        expected_arrow = "left" if ref.strand == "-" else "right"
        arrow_ok = row["arrow_direction"] == expected_arrow
        arrow_matches += int(arrow_ok)
        if not arrow_ok and len(failure_examples) < 5:
            failure_examples.append(f"arrow:{ref.qname}:{ref.chrom}:{ref.start}")

        if ref.end is not None:
            end_tested += 1
            end_ok = int(row["end"]) == ref.end
            end_matches += int(end_ok)
            if not end_ok and len(failure_examples) < 5:
                failure_examples.append(f"end:{ref.qname}:{ref.chrom}:{ref.start}")

    # PAF-to-intermediate matching and query/render checks.
    paf_map = occurrence_map(pafs, paf_key)
    row_paf_map = occurrence_map(rows, tsv_strong_key)

    query_tested = query_matches = 0
    layout_query_tested = layout_query_matches = 0
    render_seq_tested = render_seq_matches = 0
    render_events_tested = render_events_matching = 0
    symbols_tested = symbols_matching = 0
    labels_tested = labels_matching = 0

    for k, paf in paf_map.items():
        row = row_paf_map.get(k)
        if row is None:
            continue
        query_tested += 1
        q_ok = (
            int(row["query_start"]) == paf.query_start
            and int(row["query_end"]) == paf.query_end
            and int(row["query_length"]) == paf.query_length
        )
        query_matches += int(q_ok)
        # Raw query_start/query_end differences are retained internally for
        # diagnostics but are not treated as failures here. BamSnap-LRS layout
        # uses the clipping-aware layout_query_start/layout_query_end values below.

        layout_query_tested += 1
        lq_ok = (
            int(row["layout_query_start"]) == paf.query_start
            and int(row["layout_query_end"]) == paf.query_end
        )
        layout_query_matches += int(lq_ok)
        if not lq_ok and len(failure_examples) < 5:
            failure_examples.append(f"layout_query:{paf.qname}:{paf.chrom}:{paf.start}")

        (seq_ok, ev_t, ev_m, sym_t, sym_m, lab_t, lab_m) = compare_render_events(paf, row)
        if seq_ok is not None:
            render_seq_tested += 1
            render_seq_matches += int(seq_ok)
            render_events_tested += ev_t
            render_events_matching += ev_m
            symbols_tested += sym_t
            symbols_matching += sym_m
            labels_tested += lab_t
            labels_matching += lab_m
            if not seq_ok and len(failure_examples) < 5:
                failure_examples.append(f"render:{paf.qname}:{paf.chrom}:{paf.start}")

    # Group-level checks are performed only on records retained in the
    # BamSnap-LRS intermediate output. Excluded reference/PAF records must not
    # make an otherwise valid visible group fail membership/order checks.
    retained_pafs = [p for k, p in paf_map.items() if k in row_paf_map]

    ref_groups = defaultdict(list)
    paf_groups = defaultdict(list)
    row_groups = defaultdict(list)
    for r in retained_refs:
        ref_groups[r.qname].append(r)
    for p in retained_pafs:
        paf_groups[p.qname].append(p)
    for r in rows:
        row_groups[r["qname"]].append(r)

    multi_qnames = sorted(q for q, members in ref_groups.items() if len(members) > 1)
    multi_groups = len(multi_qnames)
    intermediate_multi_groups = sum(1 for members in row_groups.values() if len(members) > 1)
    group_membership_matches = 0
    group_size_matches = 0
    query_order_groups_tested = 0
    query_order_group_matches = 0
    adjacent_row_group_matches = 0
    color_groups_tested = 0
    color_group_matches = 0
    opacity_groups_tested = 0
    opacity_group_matches = 0
    body_style_groups_tested = 0
    body_style_group_matches = 0

    for q in multi_qnames:
        ref_members = ref_groups[q]
        row_members = row_groups.get(q, [])
        paf_members = paf_groups.get(q, [])

        # Membership/group-size: use the same identity strength as the reference.
        ref_member_keys = Counter(ref_key(r, strong_reference) for r in ref_members)
        row_member_keys = Counter(
            (tsv_strong_key(r) if strong_reference else tsv_weak_key(r)) for r in row_members
        )
        membership_ok = ref_member_keys == row_member_keys
        group_membership_matches += int(membership_ok)

        size_ok = (
            len(row_members) == len(ref_members)
            and all(int(r["qname_group_size"]) == len(ref_members) for r in row_members)
        )
        group_size_matches += int(size_ok)

        if len(paf_members) == len(row_members) == len(ref_members) and len(row_members) > 1:
            query_order_groups_tested += 1
            expected = sorted(paf_members, key=lambda p: (p.query_start, p.query_end, p.start, p.end))
            observed = sorted(row_members, key=lambda r: int(r["query_order_in_group"]))
            order_ok = all(
                (o["chrom"], int(o["start"]), int(o["end"]), o["strand"], int(o["mapq"]))
                == (p.chrom, p.start, p.end, p.strand, p.mapq)
                for o, p in zip(observed, expected)
            ) and [int(r["query_order_in_group"]) for r in observed] == list(range(1, len(observed) + 1))
            query_order_group_matches += int(order_ok)

            rows_in_order = [int(r["stack_row"]) for r in observed]
            adjacent_ok = rows_in_order == list(range(rows_in_order[0], rows_in_order[0] + len(rows_in_order)))
            adjacent_row_group_matches += int(adjacent_ok)

            # Same group color for all visible pieces.
            colors = [r["supplementary_group_color"] for r in observed]
            color_groups_tested += 1
            color_ok = all(c not in {"", "."} for c in colors) and len(set(colors)) == 1
            color_group_matches += int(color_ok)

            # Opacity follows renderer's query-order rule.
            expected_ops = expected_opacities(len(observed))
            observed_ops = [parse_optional_float(r["supplementary_group_opacity"]) for r in observed]
            opacity_groups_tested += 1
            opacity_ok = all(float_equal(a, b) for a, b in zip(observed_ops, expected_ops))
            opacity_group_matches += int(opacity_ok)

            # Visible body should use the same group color/opacity.
            body_style_groups_tested += 1
            body_ok = all(
                r["read_body_color"] == r["supplementary_group_color"]
                and float_equal(
                    parse_optional_float(r["read_body_opacity"]),
                    parse_optional_float(r["supplementary_group_opacity"]),
                )
                for r in observed
            )
            body_style_group_matches += int(body_ok)

            if len(failure_examples) < 5:
                if not membership_ok:
                    failure_examples.append(f"group_membership:{q}")
                elif not order_ok:
                    failure_examples.append(f"query_order:{q}")
                elif not adjacent_ok:
                    failure_examples.append(f"adjacent_rows:{q}")
                elif not color_ok:
                    failure_examples.append(f"group_color:{q}")
                elif not opacity_ok:
                    failure_examples.append(f"opacity:{q}")
                elif not body_ok:
                    failure_examples.append(f"body_style:{q}")

    # Layout collision check across all intermediate records.
    overlapping_pairs = 0
    same_row_collisions = 0
    for i in range(len(rows)):
        a = rows[i]
        for j in range(i + 1, len(rows)):
            b = rows[j]
            if a["chrom"] != b["chrom"]:
                continue
            if intervals_overlap(int(a["start"]), int(a["end"]), int(b["start"]), int(b["end"])):
                overlapping_pairs += 1
                if int(a["stack_row"]) == int(b["stack_row"]):
                    same_row_collisions += 1
                    if len(failure_examples) < 5:
                        failure_examples.append(f"collision:{a['qname']}|{b['qname']}")

    summary = {
        # provenance / input counts
        # "reference_has_cigar": strong_reference,
        "expected_records_number": len(refs),
        # "paf_alignments": len(pafs),
        "intermediate_records_number": len(rows),

        # selection / identity
        # "reference_intermediate_matches": matched_ref,
        "excluded_records": excluded_ref,
        # "unexpected_in_intermediate": unexpected_ref,
        # "paf_intermediate_matches": matched_paf,

        # BAM/SAM attributes
        # "flag_status_tested": flag_tested,
        "flag_status_matches": flag_matches,
        # "reference_end_tested": end_tested,
        # "reference_end_matches": end_matches,
        # "orientation_tested": arrow_tested,
        "orientation_matches": arrow_matches,

        # query coordinates
        # Raw query coordinates are still calculated for diagnostics, but are
        # not included in failure_examples because layout uses the clipping-aware span.
        # "raw_query_coordinates_tested": query_tested,
        # "raw_query_coordinates_matches": query_matches,
        # "query_coordinates_frompaf_tested": layout_query_tested,
        "query_coordinates_fromBamsnapLRS_matches": layout_query_matches,

        # grouping / layout on the retained subset
        "expected_multi-record_group": multi_groups,
        "intermediate_multi-record_group": intermediate_multi_groups,
        # "group_membership_matches": group_membership_matches,
        "group_size_matches": group_size_matches,
        # "query_order_tested": query_order_groups_tested,
        "query_order_matches": query_order_group_matches,
        # "adjacent_row_group_matches": adjacent_row_group_matches,
        "overlapping_record_pairs": overlapping_pairs,
        "same_row_collisions": same_row_collisions,

        # visual semantics
        # "group_color_groups_tested": color_groups_tested,
        "group_color_matches": color_group_matches,
        # "opacity_groups_tested": opacity_groups_tested,
        "opacity_order_matches": opacity_group_matches,
        # "body_style_groups_tested": body_style_groups_tested,
        # "body_style_matches": body_style_group_matches,
        # "render_event_sequences_tested": render_seq_tested,
        # "render_event_sequence_matches": render_seq_matches,
        # "render_events_tested": render_events_tested,
        # "render_events_matching": render_events_matching,
        # "symbols_tested": symbols_tested,
        # "symbols_matching": symbols_matching,
        # "labels_tested": labels_tested,
        # "labels_matching": labels_matching,

        # troubleshooting only
        "failure_examples": ";".join(failure_examples) if failure_examples else ".",
    }

    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)

    print(f"Wrote validation summary to {args.out}")
    print(
        f"selection: {matched_ref}/{len(refs)} reference alignments retained; "
        f"excluded={excluded_ref}, unexpected={unexpected_ref}"
    )
    print(
        f"flags: {flag_matches}/{flag_tested}; "
        f"layout query coordinates: {layout_query_matches}/{layout_query_tested}; "
        f"multi-groups: expected={multi_groups}, intermediate={intermediate_multi_groups}; "
        f"collisions: {same_row_collisions}/{overlapping_pairs} overlap pairs"
    )
    if failure_examples:
        print("failure examples: " + "; ".join(failure_examples))


if __name__ == "__main__":
    main()
