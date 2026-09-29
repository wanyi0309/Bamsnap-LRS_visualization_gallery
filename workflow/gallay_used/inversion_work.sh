# Dataset download
##ground truth
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/T2T-CHM13/variants_T2T-CHM13_sv_inv_sym_HGSVC2024v1.0.vcf.gz

##HiFi
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190830_220126.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190901_095311.Q20.fastq

##ONT
aria2c -x 16 -s 16 -c https://ont-open-data.s3.amazonaws.com/gm24385_2023.12/all_pass.vhg002v1.bam
samtools fastq all_pass.vhg002v1.bam > all_pass.vhg002v1.fastq

##reference genome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip -d chm13v2.0.fa.gz

#mapping
##HIFI
minimap2 -a -t 64 -x map-hifi -Q --eqx --secondary=no -K4G --MD -Y -L chm13v2.0.fa.gz m64011_190830_220126.Q20.fastq m64011_190901_095311.Q20.fastq | samtools view -bh -@ 64 - | samtools sort -@ 64 -o HIFI.primary.sort.bam -
samtools index HIFI.primary.sort.bam
##ONT
minimap2 -a -t 64 -x map-ont -Q --eqx --secondary=no -K4G --MD -Y -L chm13v2.0.fa.gz all_pass.vhg002v1.fastq | samtools view -bh -@ 64 - | samtools sort -@ 64 -o ONT.primary.sort.bam -
samtools index ONT.primary.sort.bam

# Variant selection
bcftools view -s NA24385 -Ou variants_T2T-CHM13_sv_inv_sym_HGSVC2024v1.0.vcf.gz | bcftools view -i 'GT="alt" && INFO/SVLEN<100000' | grep -v "#" > inversion.vcf

#svlen
awk 'BEGIN{OFS="\t"} {svlen=0
    n=split($8,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^SVLEN=/){
            split(a[i],b,"=")
            svlen=b[2]+0
            if(svlen<0) svlen=-svlen}}
    print $1,$2-1,$2+svlen-1,"INV"NR}' inversion.vcf  > inversion.padding0.bed
    
awk 'BEGIN{OFS="\t"} {svlen=0
    n=split($8,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^SVLEN=/){
            split(a[i],b,"=")
            svlen=b[2]+0
            if(svlen<0) svlen=-svlen}}
    print $1,$2-1-2000,$2+svlen-1+2000,"INV"NR}' inversion.vcf  > inversion.padding2000.bed


#Visualization
time -v bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions inversion.vcf --out-prefix bamsnap/inversion.png --fa chm13v2.0.fa  --show-axis --show-coverage --padding 2000  --show-supp > bamsnap_inv.runtime.log 2>&1 

time -v hawkeye.py sv_browse -g hg_38 -i HIFI.primary.sort.bam,ONT.primary.sort.bam  -b inversion.padding0.bed -r chm13v2.0.fa -o svhawkeye/ -q 0 -I 0 -d 2000 -F png -f bed --sv_min_length=0  > svhawkeye_inv.runtime.log 2>&1 

time -v wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -u -g chm13v2.0.fa -R inversion.padding2000.bed HIFI.primary.sort.bam ONT.primary.sort.bam  > wally_inv.runtime.log 2>&1 

tar  -czvf  wally.INV.tar.gz  INV*png
tar  -czvf bamsnap_lrs.INV.tar.gz bamsnap/*png
tar  -czvf svhawkeye.INV.tar.gz svhawkeye/figure/*png
