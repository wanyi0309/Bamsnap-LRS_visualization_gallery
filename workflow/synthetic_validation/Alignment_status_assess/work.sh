## BAM source
HIFI.primary.sort.bam ONT.primary.sort.bam from https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/inversion_work.sh
SRR8955953.pacbio.hg38.sorted.bam SRR16005301.ont.hg38.sorted.bam from https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/duplication_work.sh


cat INV.check_region.txt | while read region CASE
do
    mkdir -p $CASE && cd $CASE
    samtools view -h HIFI.primary.sort.bam $region | k8 paftools.js sam2paf - > $CASE.HIFI.paf
    samtools view HIFI.primary.sort.bam $region | awk '{print $1"\t"$2 "\t"$3 "\t"$4 "\t"$5 }' > $CASE.HIFI.sam_reference.tsv
    python dump_bamsnap_intermediate.py --bamsnap-src path/to/bamsnap-lrs/src/ --bam HIFI.primary.sort.bam  --pos $region --show-supp  --show-secondary  --detail high --overview-detail show --out $CASE.HIFI.bamsnap_intermediate.tsv
    python validate_bamsnap_intermediate.py --reference $CASE.HIFI.sam_reference.tsv --paf $CASE.HIFI.paf --intermediate $CASE.HIFI.bamsnap_intermediate.tsv -o $CASE.HIFI.validation.tsv
    cd ..
done
cat INV.check_region.txt | while read region CASE
do
    mkdir -p $CASE && cd $CASE
    samtools view -h ONT.primary.sort.bam $region | k8 paftools.js sam2paf - > $CASE.ONT.paf
    samtools view ONT.primary.sort.bam $region | awk '{print $1"\t"$2 "\t"$3 "\t"$4 "\t"$5 }' > $CASE.ONT.sam_reference.tsv
    python dump_bamsnap_intermediate.py --bamsnap-src path/to/bamsnap-lrs/src/ --bam ONT.primary.sort.bam  --pos $region --show-supp  --show-secondary  --detail high --overview-detail show --out $CASE.ONT.bamsnap_intermediate.tsv
    python validate_bamsnap_intermediate.py --reference $CASE.ONT.sam_reference.tsv --paf $CASE.ONT.paf --intermediate $CASE.ONT.bamsnap_intermediate.tsv -o $CASE.ONT.validation.tsv
    cd ..
done

cat DUP.check_region.txt | while read region CASE
do
    mkdir -p $CASE && cd $CASE
    samtools view -h SRR8955953.pacbio.hg38.sorted.bam $region | k8 paftools.js sam2paf - > $CASE.HIFI.paf
    samtools view SRR8955953.pacbio.hg38.sorted.bam $region | awk '{print $1"\t"$2 "\t"$3 "\t"$4 "\t"$5 }' > $CASE.HIFI.sam_reference.tsv
    python dump_bamsnap_intermediate.py --bamsnap-src path/to/bamsnap-lrs/src/ --bam SRR8955953.pacbio.hg38.sorted.bam  --pos $region --show-supp --show-secondary --detail high --overview-detail show --out $CASE.HIFI.bamsnap_intermediate.tsv
    python validate_bamsnap_intermediate.py --reference $CASE.HIFI.sam_reference.tsv --paf $CASE.HIFI.paf --intermediate $CASE.HIFI.bamsnap_intermediate.tsv -o $CASE.HIFI.validation.tsv
    cd ..
done
cat DUP.check_region.txt | while read region CASE
do
    mkdir -p $CASE && cd $CASE
    samtools view -h SRR16005301.ont.hg38.sorted.bam $region | k8 paftools.js sam2paf - > $CASE.ONT.paf
    samtools view SRR16005301.ont.hg38.sorted.bam $region | awk '{print $1"\t"$2 "\t"$3 "\t"$4 "\t"$5 }' > $CASE.ONT.sam_reference.tsv
    python dump_bamsnap_intermediate.py --bamsnap-src path/to/bamsnap-lrs/src/ --bam SRR16005301.ont.hg38.sorted.bam  --pos $region --show-supp --show-secondary --detail high --overview-detail show --out $CASE.ONT.bamsnap_intermediate.tsv
    python validate_bamsnap_intermediate.py --reference $CASE.ONT.sam_reference.tsv --paf $CASE.ONT.paf --intermediate $CASE.ONT.bamsnap_intermediate.tsv -o $CASE.ONT.validation.tsv
    cd ..
done
