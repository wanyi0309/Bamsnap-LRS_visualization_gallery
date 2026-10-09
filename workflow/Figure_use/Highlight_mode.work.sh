# For details on downloading and generating the BAM, FASTA and VCF files used in this script
# please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/phase_work.sh

# select region: chr12	74994312	75006232	random_region6

bamsnap-lrs highlight --highlight-vcf HG002.phased_het_snps.vcf.gz --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam --pos chr12:74994312-75006232  --out bamsnap_SNV.svg --fa hg38.fa  --no-hap-filter --padding 0 --mapq 0 

wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -g hg38.fa -r random_50_regions.forwally.bed HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam

samtools view HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam chr12:74994312-75006232 -bh > highlight.bam
samtools index highlight.bam
