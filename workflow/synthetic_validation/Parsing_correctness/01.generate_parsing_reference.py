#!/usr/bin/env python3
import argparse
import random
from pathlib import Path

def write_fasta(path: Path, name: str, seq: str, line_width: int = 60):
    with path.open("w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i:i + line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Generate the synthetic reference used for BamSnap-LRS parsing validation.")
    parser.add_argument("-o", "--output", default="synthetic_reference.fa", help="Output FASTA path (default: synthetic_reference.fa)")
    parser.add_argument("--length", type=int, default=1000, help="Reference length in bp (default: 1000)")
    parser.add_argument("--seed", type=int, default=20260910, help="Random seed (default: 20260910)")
    parser.add_argument("--chrom_name", default="chrTest", help="Synthetic reference chromosome name")
    parser.add_argument("--dup_from", default=None, help="Source duplicated interval (format: start-end), 0-based half-open (optional)")
    parser.add_argument("--dup_to", default=None, help="Destination duplicated interval (format: start-end), 0-based half-open (optional)")
    args = parser.parse_args()

    do_dup = args.dup_from is not None or args.dup_to is not None
    if do_dup:
        if args.dup_from is None or args.dup_to is None:
            raise ValueError("--dup_from and --dup_to must be provided together, or neither.")

        try:
            dup_from_start, dup_from_end = map(int, args.dup_from.split("-"))
            dup_to_start, dup_to_end = map(int, args.dup_to.split("-"))
        except ValueError:
            raise ValueError("--dup_from and --dup_to must use START-END format, e.g. 100-160")

        if dup_from_start < 0 or dup_to_start < 0:
            raise ValueError("Duplicate interval coordinates must be non-negative.")
        if dup_from_end <= dup_from_start or dup_to_end <= dup_to_start:
            raise ValueError("Duplicate interval end must be greater than start.")
        if (dup_from_end - dup_from_start) != (dup_to_end - dup_to_start):
            raise ValueError("Source and destination duplicate intervals must have the same length.")
        if dup_from_end > args.length or dup_to_end > args.length:
            raise ValueError(f"Duplicate intervals must fall within the reference length ({args.length} bp).")

    rng = random.Random(args.seed)
    seq = [rng.choice("ACGT") for _ in range(args.length)]

    if do_dup:
        seq[dup_to_start:dup_to_end] = seq[dup_from_start:dup_from_end]
    seq = "".join(seq)

    output = Path(args.output)
    write_fasta(output, args.chrom_name, seq)

    assert len(seq) == args.length
    if do_dup:
        assert seq[dup_from_start:dup_from_end] == seq[dup_to_start:dup_to_end]

    print(f"Generated: {output}")
    print(f"Reference name: {args.chrom_name}")
    print(f"Reference length: {len(seq)} bp")
    print(f"Random seed: {args.seed}")
    if do_dup:
        print(f"Duplicated segment: {args.chrom_name}[{dup_from_start}:{dup_from_end}] == {args.chrom_name}[{dup_to_start}:{dup_to_end}] (0-based, half-open)")
        print(f"Equivalent samtools coordinates: {args.chrom_name}:{dup_from_start + 1}-{dup_from_end} == {args.chrom_name}:{dup_to_start + 1}-{dup_to_end}")
    else:
        print("Duplicated segment: none (no --dup_from/--dup_to provided)")
    print(f"Next step: samtools faidx {output}")

if __name__ == "__main__":
    main()