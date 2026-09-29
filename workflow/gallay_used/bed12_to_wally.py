#!/usr/bin/env python3
import argparse
import gzip
import sys


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.lower().endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def chrom_key(chrom):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return (0, int(c), "")
    rank = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if c in rank:
        return (0, rank[c], "")
    return (1, 0, chrom)


def parse_bed12(line, line_no):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        raise ValueError(f"line {line_no}: expected >=12 BED columns")

    chrom = f[0]
    start = int(f[1])
    end = int(f[2])
    name = f[3] or "."
    strand = f[5] if f[5] in {"+", "-"} else "."
    block_count = int(f[9])
    sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
    starts = [int(x) for x in f[11].rstrip(",").split(",") if x]

    if len(sizes) != block_count or len(starts) != block_count:
        raise ValueError(
            f"line {line_no}: blockCount={block_count}, "
            f"sizes={len(sizes)}, starts={len(starts)}"
        )

    exons = [(start + rel, start + rel + size)
             for rel, size in zip(starts, sizes)]
    return chrom, start, end, name, strand, block_count, exons


def is_3bp_marker(start, end, block_count):
    return block_count == 1 and (end - start) == 3


def main():
    ap = argparse.ArgumentParser(
        description="Convert BED12/BED12+ transcript models to Wally-style BED6."
    )
    ap.add_argument("input", help="Input BED12/BED12+ (.bed/.bed.gz/.bgz)")
    ap.add_argument("-o", "--output", default="-", help="Output BED6 [stdout]")
    ap.add_argument("--keep-3bp-markers", action="store_true",
                    help="Keep blockCount=1, span=3 bp records (excluded by default)")
    args = ap.parse_args()

    rows = []
    total = 0
    excluded = 0
    retained = 0

    with open_text(args.input) as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue

            chrom, tx_start, tx_end, name, strand, block_count, exons = \
                parse_bed12(line, line_no)
            total += 1

            if is_3bp_marker(tx_start, tx_end, block_count) and not args.keep_3bp_markers:
                excluded += 1
                continue

            retained += 1
            n = len(exons)

            for genomic_idx, (exon_start, exon_end) in enumerate(exons):
                exon_no = n - genomic_idx if strand == "-" else genomic_idx + 1
                rows.append((
                    chrom, exon_start, exon_end,
                    f"{name} ({exon_no}/{n})", "exon", strand, 0
                ))

            rows.append((
                chrom, tx_start, tx_end, name, "transcript", strand, 1
            ))

    rows.sort(key=lambda r: (
        chrom_key(r[0]), r[1], r[2], r[6], r[3]
    ))

    out = sys.stdout if args.output == "-" else open(args.output, "wt")
    try:
        for chrom, start, end, label, feature_type, strand, _ in rows:
            out.write(
                f"{chrom}\t{start}\t{end}\t{label}\t{feature_type}\t{strand}\n"
            )
    finally:
        if out is not sys.stdout:
            out.close()

    print(f"Input BED12 records: {total}", file=sys.stderr)
    print(f"Excluded 3-bp single-block records: {excluded}", file=sys.stderr)
    print(f"Retained transcript models: {retained}", file=sys.stderr)


if __name__ == "__main__":
    main()
