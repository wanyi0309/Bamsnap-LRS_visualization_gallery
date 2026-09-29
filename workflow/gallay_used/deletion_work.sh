#Dataset download
##ground truth
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed

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
bcftools view -H -R HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed   -i 'INFO/SVTYPE="DEL" && INFO/SVLEN>50' HG002_CHM13v2.0_v5.0q_stvar.vcf.gz | shuf --random-source=<(yes 123456789) -n 50 > deletion_GIAB_random50.vcf

awk 'BEGIN{OFS="\t"} {svlen=0
    n=split($8,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^SVLEN=/){
            split(a[i],b,"=")
            svlen=b[2]+0
            if(svlen<0) svlen=-svlen}}
    print $1,$2-1,$2+svlen-1,"DEL"NR}' deletion_GIAB_random50.vcf  > deletion_GIAB_random50.bed

awk 'BEGIN{OFS="\t"} {svlen=0
    n=split($8,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^SVLEN=/){
            split(a[i],b,"=")
            svlen=b[2]+0
            if(svlen<0) svlen=-svlen}}
    print $1,$2-1-500,$2+svlen-1+500,"DEL"NR}' deletion_GIAB_random50.vcf  > deletion_GIAB_random50.padding500.bed

#Visualization
time -v python bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions deletion_GIAB_random50.vcf --out-prefix bamsnap/deletion_ONT_HIFI.png --fa chm13v2.0.fa  --show-axis --show-coverage --padding 500  --detail high --overview-detail show > bamsnap_del.runtime.log 2>&1 

time -v hawkeye.py sv_browse -g chm13 -i HIFI.primary.sort.bam,ONT.primary.sort.bam  -b deletion_GIAB_random50.padding500.bed -r chm13v2.0.fa -o svhawkeye/ -q 0 -I 0 -d 500 -F png -f bed --sv_min_length=0  > hawkeye_del.runtime.log 2>&1 

time -v wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -g chm13v2.0.fa -R deletion_GIAB_random50.padding500.bed HIFI.primary.sort.bam ONT.primary.sort.bam > wally_del.runtime.log 2>&1 

tar  -czvf  wally.DEL.tar.gz  DEL*png
tar  -czvf bamsnap_lrs.DEL.tar.gz bamsnap/*png
tar  -czvf svhawkeye.DEL.tar.gz svhawkeye/figure/*png
