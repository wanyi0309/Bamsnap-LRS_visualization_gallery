# For details on downloading and generating the BAM and FASTA files used in this script
# please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/duplication_work.sh

## select region: chr10	74808908	74818815	DUP13

bamsnap-lrs dna --bam SRR8955953.pacbio.hg38.sorted.bam --pos chr10:74808908-74818815 --out bamsnap_duplication.svg --fa hg38.fa  --show-axis --show-coverage --padding 0 --show-supp  --width 700 

echo -e "chr10\t74808908\t74818815\tDUP13" > target.bed 
hawkeye.py sv_browse -g hg_38 -i SRR8955953.pacbio.hg38.sorted.bam  -b target.bed -r hg38.fa -o svhawkeye/ -q 0 -I 0 -d 0 -F pdf -f bed --sv_min_length=0  

wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 700 -u -g hg38.fa -r chr10:74808908-74818815 SRR8955953.pacbio.hg38.sorted.bam 
 
samtools view SRR8955953.pacbio.hg38.sorted.bam chr10:74808908-74818815 -bh > dup.bam
samtools index dup.bam
