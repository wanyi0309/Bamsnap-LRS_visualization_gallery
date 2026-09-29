#Dataset download
##ground truth
wget https://cgl.gi.ucsc.edu/data/LRGASP/annotations/human/LRGASP_manual_annotation.human.cDNA_PacBio.bed.gz

##BAM
wget https://www.encodeproject.org/files/ENCFF985LGZ/@@download/ENCFF985LGZ.bam
wget https://www.encodeproject.org/files/ENCFF373TKM/@@download/ENCFF373TKM.bam
wget https://www.encodeproject.org/files/ENCFF388HXU/@@download/ENCFF388HXU.bam

##reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna.gz hg38.fa.gz
gzip -d hg38.fa.gz

#Data processing
# https://cgl.gi.ucsc.edu/data/LRGASP/annotations/human/LRGASP_manual_annotation.human.cDNA_PacBio.bed.gz
gzip -dc LRGASP_manual_annotation.human.cDNA_PacBio.bed.gz | cut -f1-12   > LRGASP_manual_annotation.human.cDNA_PacBio.bed
gzip -dc LRGASP_manual_annotation.human.cDNA_PacBio.bed.gz | awk '$11!="3,"' | cut -f1-12 > LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed
python bed12_to_wally.py LRGASP_manual_annotation.human.cDNA_PacBio.bed -o LRGASP_manual_annotation.human.cDNA_PacBio.wally.bed
sort -k1,1V -k2,2n -k3,3n LRGASP_manual_annotation.human.cDNA_PacBio.wally.bed > LRGASP_manual_annotation.human.cDNA_PacBio.wally.sorted.bed
bgzip LRGASP_manual_annotation.human.cDNA_PacBio.wally.sorted.bed
tabix -p bed  LRGASP_manual_annotation.human.cDNA_PacBio.wally.sorted.bed.gz

# https://www.encodeproject.org/experiments/ENCSR507JOF/ human WTC11 hg38
/share/app/samtools/1.11/bin/samtools sort -@ 16 -o ENCFF985LGZ.sorted.bam ENCFF985LGZ.bam
/share/app/samtools/1.11/bin/samtools sort -@ 16 -o ENCFF373TKM.sorted.bam ENCFF373TKM.bam
/share/app/samtools/1.11/bin/samtools sort -@ 16 -o ENCFF388HXU.sorted.bam ENCFF388HXU.bam
/share/app/samtools/1.11/bin/samtools merge -@ 16 human.hg38.Pacbio.merged.bam ENCFF985LGZ.sorted.bam ENCFF373TKM.sorted.bam ENCFF388HXU.sorted.bam
/share/app/samtools/1.11/bin/samtools index human.hg38.Pacbio.merged.bam

## plot region
python bed12_gene_regions.py LRGASP_manual_annotation.human.cDNA_PacBio.bed.gz --gene-map Supplementary_Data14_human_genes.tsv --mapped-only -o LRGASP_cDNA_PacBio_50genes_regions.bed --details LRGASP_cDNA_PacBio_50genes_regions.tsv
awk '{print $1"\t"($2-1000)"\t"($3+1000)"\t"$4}' LRGASP_cDNA_PacBio_50genes_regions.bed > LRGASP_cDNA_PacBio_50genes_regions.padding.bed

#Visualization
time -v bamsnap-lrs rna --bam human.hg38.Pacbio.merged.bam --regions LRGASP_cDNA_PacBio_50genes_regions.bed --out-prefix bamsnap/RNA.png  --mapq 0 --fa hg38.fa  --padding 1000  --bed LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed --show-axis > bamsnap_RNA.runtime.log  2>&1 

time -v hawkeye.py rna_browse -g hg38 -i human.hg38.Pacbio.merged.bam  -b LRGASP_cDNA_PacBio_50genes_regions.padding.bed -o svhawkeye_test -q 0 -I 0  -F png > svhawkeye_RNA.runtime.log  2>&1 

time -v wally-v0.9.2-linux-amd64 region -g hg38.fa -R LRGASP_cDNA_PacBio_50genes_regions.padding.bed human.hg38.Pacbio.merged.bam --map-qual 0 --bed LRGASP_manual_annotation.human.cDNA_PacBio.wally.sorted.bed.gz > wally_RNA.runtime.log   2>&1 

tar  -czvf  wally.RNA.tar.gz  *png
tar  -czvf bamsnap_lrs.RNA.tar.gz bamsnap/*png
tar  -czvf svhawkeye.RNA.tar.gz svhawkeye_test/figure/*png
