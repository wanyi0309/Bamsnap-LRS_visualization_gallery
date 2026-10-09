# Dataset download
##ground truth
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz

##HiFi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam

##reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna.gz hg38.fa.gz
gzip -d hg38.fa.gz
samtools faidx hg38.fa


# region processing
awk 'BEGIN{srand(12345);for(i=1;i<=50;i++){len=int(1000 + rand()*(50000-1000+1)); print "chr1\t0\t"len } }' > random_length_seed.bed
grep '^chr' hg38.fa.fai | cut -f1-2 > hg38.genome
# hg38.chromosome_band.txt down from https://genome.ucsc.edu/cgi-bin/hgTables
awk '$5=="gvar" || $5=="acen" || $5=="stalk"' hg38.chromosome_band.txt > hg38.chromosome_band.gvar_acen_stalk.txt

/share/app/bedtools/2.29.2/bin/bedtools shuffle -seed 123456789 -i random_length_seed.bed -g hg38.genome -excl hg38.chromosome_band.gvar_acen_stalk.txt > random_50_regions.bed
awk '{print $1"\t"$2"\t"$3"\trandom_region"NR}' random_50_regions.bed > random_50_regions.forwally.bed


# Variant selection
/share/app/bcftools/1.11/bin/bcftools view -p -g het -i 'TYPE="snp" && GT!~"\."' HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz -Oz -o HG002.phased_het_snps.vcf.gz
/share/app/bcftools/1.11/bin/bcftools index HG002.phased_het_snps.vcf.gz

# Plot
time -v bamsnap-lrs highlight --highlight-vcf HG002.phased_het_snps.vcf.gz --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam --regions random_50_regions.bed  --out-prefix bamsnap/SNV.png --fa hg38.fa  --no-hap-filter --padding 0 --mapq 0 > bamsnap_hap.log 2>&1

time -v wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -ghg38.fa -R random_50_regions.forwally.bed HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam > wally_hap.log 2>&1

tar  -czvf  wally.Highlight.tar.gz  random_region*png
tar  -czvf bamsnap_lrs.Highlight.tar.gz bamsnap/*png
