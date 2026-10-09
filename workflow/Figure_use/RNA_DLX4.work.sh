# For details on downloading and generating the BAM, FASTA, and annotation files used in this script
# please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/RNA_work.sh

## select region: chr17	49967564	49975971	DLX4

bamsnap-lrs rna --bam human.hg38.Pacbio.merged.bam --pos chr17:49967564-49975971 --out bamsnap_RNA.svg  --mapq 0 --fa hg38.fa  --padding 1000  --bed LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed --show-axis  --width 800 --read-height 12

echo -e "chr17\t49966564\t49976971\tDLX4" > target.bed # region padding 1000
hawkeye.py rna_browse -g hg38 -i human.hg38.Pacbio.merged.bam  -b target.bed -o svhawkeye_test -q 0 -I 0  -F pdf 

wally-v0.9.2-linux-amd64 region -g hg38.fa -R target.bed human.hg38.Pacbio.merged.bam --map-qual 0 --bed LRGASP_manual_annotation.human.cDNA_PacBio.wally.sorted.bed.gz  -x 800


samtools view human.hg38.Pacbio.merged.bam chr17:49966564-49976971 -bh > RNA.bam
samtools index RNA.bam
