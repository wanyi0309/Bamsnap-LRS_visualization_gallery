#Dataset download
##ground truth
# we obtained duplications from a SV high confidence consensus call set from Talsania et al. 2022 [10.1186/s13059-022-02816-6] 
# and sample 50 cases (select DUPs with a length of less than 50 kb from the table https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13059-022-02816-6/MediaObjects/13059_2022_2816_MOESM4_ESM.xlsx)
# For details on the selected DUPs, please refer to duplication.bed

##pacbio
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8955953/SRR8955953
fastq-dump SRR8955953 --split-3 --gzip --defline-qual '+'  -O fq
##ONT
wget 	https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR16005301/SRR16005301
fastq-dump SRR16005301 --split-3 --gzip --defline-qual '+'  -O fq

##reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna.gz hg38.fa.gz
gzip -d hg38.fa.gz

#Mapping
##ONT
mkdir -p /data/work/00.data/SRR16005301/bam
minimap2 -ax map-ont --MD -t 16 hg38.fa SRR16005301.fastq.gz | samtools view -@ 16 -bS - | samtools sort -m 16G -@ 16 -o /data/work/00.data/SRR16005301/bam/SRR16005301.ont.hg38.sorted.bam -
samtools index -@ 16 /data/work/00.data/SRR16005301/bam/SRR16005301.ont.hg38.sorted.bam
##pacbio
mkdir -p /data/work/00.data/SRR8955953/bam
minimap2 -ax map-pb --MD -t 16 hg38.fa /data/work/00.data/SRR8955953/SRR8955953.fastq.gz | samtools view -@ 8 -bS - | samtools sort -m 4G -@ 8 -o /data/work/00.data/SRR8955953/bam/SRR8955953.pacbio.hg38.sorted.bam -
samtools index -@ 8 /data/work/00.data/SRR8955953/bam/SRR8955953.pacbio.hg38.sorted.bam

awk '{print $1"\t"($2-2000)"\t"($3+2000)"\t"$4NR}' duplication.bed > duplication.padding2000.bed

#Visualization
time -v python /data/work/01.bamsnap_lrs/12.modify.v11.0920/bin/bamsnap-lrs dna --bam /data/work/00.data/SRR16005301/bam/SRR16005301.ont.hg38.sorted.bam /data/work/00.data/SRR8955953/bam/SRR8955953.pacbio.hg38.sorted.bam --regions duplication.bed --out-prefix bamsnap/duplication_ONT.png --fa hg38.fa  --show-axis --show-coverage --padding 2000  --show-supp  > bamsnap_dup.runtime.log 2>&1 

time -v hawkeye.py sv_browse -g hg_38 -i /data/work/00.data/SRR16005301/bam/SRR16005301.ont.hg38.sorted.bam,/data/work/00.data/SRR8955953/bam/SRR8955953.pacbio.hg38.sorted.bam  -b duplication.bed -r hg38.fa -o svhawkeye/ -q 0 -I 0 -d 2000 -F png -f bed --sv_min_length=0  > svhawkeye_dup.runtime.log 2>&1 

time -v wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -u -g hg38.fa -R duplication.padding2000.bed /data/work/00.data/SRR16005301/bam/SRR16005301.ont.hg38.sorted.bam /data/work/00.data/SRR8955953/bam/SRR8955953.pacbio.hg38.sorted.bam > wally_dup.runtime.log 2>&1 
 
tar  -czvf  wally.DUP.tar.gz  DUP*png
tar  -czvf bamsnap_lrs.DUP.tar.gz bamsnap/*png
tar  -czvf svhawkeye.INV.DUP.gz svhawkeye/figure/*png
