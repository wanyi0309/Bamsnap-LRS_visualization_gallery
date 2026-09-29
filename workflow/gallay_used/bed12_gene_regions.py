#!/usr/bin/env python3
import argparse
import gzip
import sys
from collections import defaultdict


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.lower().endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_gene_map(path):
    """Read a two-column TSV: gene_symbol<TAB>gene_id."""
    if not path:
        return {}
    mapping = {}
    with open_text(path) as f:
        first = True
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}: line {line_no}: expected >=2 tab-separated columns")
            a, b = parts[0].strip(), parts[1].strip()
            if first:
                first = False
                if a.lower() in {"gene_symbol", "symbol", "gene"} and b.lower() in {
                    "gene_id", "geneid", "ensembl_gene_id"
                }:
                    continue
            mapping[b] = a
    return mapping


def natural_chrom_key(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return (0, int(c), "")
    rank = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if c in rank:
        return (0, rank[c], "")
    return (1, 0, chrom)


def parse_header(line):
    return {name: i for i, name in enumerate(line.rstrip("\n").lstrip("#").split("\t"))}


def is_3bp_marker(start, end, block_count):
    return block_count == 1 and (end - start) == 3


def main():
    ap = argparse.ArgumentParser(
        description="Generate the widest retained BED12/BED12+ interval per gene."
    )
    ap.add_argument("bed", help="Input BED12/BED12+ (.bed/.bed.gz/.bgz)")
    ap.add_argument("-o", "--output", required=True,
                    help="Output BED4: chrom, start, end, gene label")
    ap.add_argument("--gene-map",
                    help="Optional TSV: gene_symbol<TAB>gene_id")
    ap.add_argument("--gene-id-column", default="geneName2",
                    help="BED12+ column containing gene ID [geneName2]")
    ap.add_argument("--mapped-only", action="store_true",
                    help="Keep only gene IDs present in --gene-map")
    ap.add_argument("--details",
                    help="Optional detailed TSV")
    ap.add_argument("--keep-3bp-markers", action="store_true",
                    help="Keep blockCount=1, span=3 bp records (excluded by default)")
    args = ap.parse_args()

    gene_map = read_gene_map(args.gene_map)
    header = None
    data = defaultdict(lambda: {
        "chrom": None, "strand": None, "start": None, "end": None,
        "transcript_models": 0, "multi_exon_models": 0,
        "max_exon_count": 0, "excluded_3bp_markers": 0
    })
    total_input = 0
    total_excluded = 0

    with open_text(args.bed) as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.lower().startswith("#chrom"):
                    header = parse_header(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{args.bed}: line {line_no}: expected BED12+, got {len(fields)} columns")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5] if len(fields) > 5 else "."
            block_count = int(fields[9])

            if header is not None and args.gene_id_column in header:
                gene_idx = header[args.gene_id_column]
            elif args.gene_id_column == "geneName2" and len(fields) > 18:
                gene_idx = 18
            else:
                raise ValueError(
                    f"Cannot locate gene-ID column '{args.gene_id_column}'."
                )

            gene_id = fields[gene_idx].strip()
            if not gene_id:
                continue
            if args.mapped_only and gene_map and gene_id not in gene_map:
                continue

            total_input += 1
            rec = data[gene_id]

            if is_3bp_marker(start, end, block_count) and not args.keep_3bp_markers:
                rec["excluded_3bp_markers"] += 1
                total_excluded += 1
                continue

            if rec["chrom"] is None:
                rec["chrom"] = chrom
            elif rec["chrom"] != chrom:
                raise ValueError(f"{gene_id} appears on multiple chromosomes")

            if rec["strand"] is None:
                rec["strand"] = strand
            elif rec["strand"] in {"+", "-"} and strand in {"+", "-"} and rec["strand"] != strand:
                raise ValueError(f"{gene_id} appears on multiple strands")

            rec["start"] = start if rec["start"] is None else min(rec["start"], start)
            rec["end"] = end if rec["end"] is None else max(rec["end"], end)
            rec["transcript_models"] += 1
            rec["multi_exon_models"] += int(block_count >= 2)
            rec["max_exon_count"] = max(rec["max_exon_count"], block_count)

    rows = []
    for gene_id, rec in data.items():
        if rec["start"] is None:
            continue
        rows.append({
            "gene_symbol": gene_map.get(gene_id, gene_id),
            "gene_id": gene_id,
            "chrom": rec["chrom"],
            "start": rec["start"],
            "end": rec["end"],
            "strand": rec["strand"] or ".",
            "span_bp": rec["end"] - rec["start"],
            "transcript_models": rec["transcript_models"],
            "multi_exon_models": rec["multi_exon_models"],
            "max_exon_count": rec["max_exon_count"],
            "excluded_3bp_markers": rec["excluded_3bp_markers"],
        })

    rows.sort(key=lambda r: (
        natural_chrom_key(r["chrom"]), r["start"], r["end"], r["gene_symbol"]
    ))

    with open(args.output, "wt") as out:
        for r in rows:
            out.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['gene_symbol']}\n")

    if args.details:
        cols = [
            "gene_symbol", "gene_id", "chrom", "start", "end", "strand",
            "span_bp", "transcript_models", "multi_exon_models",
            "max_exon_count", "excluded_3bp_markers"
        ]
        with open(args.details, "wt") as out:
            out.write("\t".join(cols) + "\n")
            for r in rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"Input annotation records considered: {total_input}", file=sys.stderr)
    print(f"Excluded 3-bp single-block records: {total_excluded}", file=sys.stderr)
    print(f"Genes written: {len(rows)}", file=sys.stderr)

    if gene_map:
        observed = set(data)
        expected = set(gene_map)
        print(f"Gene-map IDs matched: {len(observed & expected)}/{len(expected)}", file=sys.stderr)


if __name__ == "__main__":
    main()
