# final select region  chr5:705000-710000
# Mc Cartney AM, Shafin K, Alonge M, Bzikadze AV, Formenti G, Fungtammasan A, Howe K, Jain C, Koren S, Logsdon GA et al: Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods 2022, 19(6):687-695.

#Dataset download
##HiFi
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190830_220126.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190901_095311.Q20.fastq

##reference genome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip -d chm13v2.0.fa.gz

#mapping
##HIFI
minimap2 -a -t 64 -x map-hifi -Q --eqx --secondary=no -K4G --MD -Y -L chm13v2.0.fa.gz m64011_190830_220126.Q20.fastq m64011_190901_095311.Q20.fastq | samtools view -bh -@ 64 - | samtools sort -@ 64 -o HIFI.primary.sort.bam -
samtools index HIFI.primary.sort.bam

bamsnap-lrs dna --bam HIFI.primary.sort.bam  --pos chr5:705000-710000 --out assembly_exam.svg   --show-coverage  --fa chm13v2.0.fa  --width 700  --show-axis --overview-detail show --mapq 0

echo -e "chr5\t705000\t710000" > target.bed
hawkeye.py sv_browse -i  HIFI.primary.sort.bam -b target.bed  -f bed -t 1 -r chm13v2.0.fa  -g chm13 -o assem_valid -d 0 -F pdf -q 0 

wally-v0.9.2-linux-amd64 region  -g chm13v2.0.fa -r chr5:705000-710000 HIFI.primary.sort.bam  -x 700 -q 0


samtools view HIFI.primary.sort.bam chr5:705000-710000 -bh > assembly_validation.bam
samtools index assembly_validation.bam
