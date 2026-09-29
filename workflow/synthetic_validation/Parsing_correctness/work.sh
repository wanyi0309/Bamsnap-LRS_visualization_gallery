## first step: generate reference genome sequence
python 01.generate_parsing_reference.py -o synthetic_reference.fa --length 1000 --seed 20260910 --chrom_name chrTest --dup_from 100-160 --dup_to 700-760

## second step:
# design custom alignment records, including file: parsing_manual_cigar_md_cs.tsv  parsing_expected_records.tsv parsing_expected_segments.tsv 

## third step:
python 02.generate_parsing_bam_manual.py -r synthetic_reference.fa -d parsing_manual_cigar_md_cs.tsv -o synthetic_parsing.bam
# /share/app/samtools/1.11/bin/samtools calmd synthetic_parsing.bam synthetic_reference.fa > synthetic_parsing.calmd.sam # double check MD 

## forth step:
python 03.test_parsing.py \
  --bam synthetic_parsing.bam \
  --fa synthetic_reference.fa \
  --records-truth parsing_expected_records.tsv \
  --segments-truth parsing_expected_segments.tsv \
  --bamsnap-src path/to/bamsnap-lrs/src \
  -o parsing_validation_results.tsv
