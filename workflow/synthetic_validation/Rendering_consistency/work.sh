## All relevant BAM files can be obtained via the download links and methods specified in `*_work.sh` located in `https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/tree/main/workflow/gallay_used`.

## insertion 
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions insertion_select.region.bed --out-prefix ins/insertion_ONT_HIFI.png --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions insertion_select.region.bed --out-prefix ins/insertion_ONT_HIFI.svg --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions insertion_select.region.bed --out-prefix ins/insertion_ONT_HIFI.pdf --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show

## deletion
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions deletion_select.region.bed --out-prefix del/deletion_ONT_HIFI.png --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions deletion_select.region.bed --out-prefix del/deletion_ONT_HIFI.svg --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions deletion_select.region.bed --out-prefix del/deletion_ONT_HIFI.pdf --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --detail high --overview-detail show

# inversion 
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions inversion_select.region.bed --out-prefix inv/inversion.png --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --show-supp
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions inversion_select.region.bed --out-prefix inv/inversion.svg --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --show-supp
bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions inversion_select.region.bed --out-prefix inv/inversion.pdf --fa /data/work/00.data/chm13v2.0.fa  --show-axis --show-coverage --padding 0  --show-supp

# duplication
bamsnap-lrs dna --bam SRR16005301.ont.hg38.sorted.bam SRR8955953.pacbio.hg38.sorted.bam --regions duplication_select.region.bed --out-prefix dup/duplication.png --fa /Files/chenjin/WGS_hg38.fa  --show-axis --show-coverage --padding 0  --show-supp
bamsnap-lrs dna --bam SRR16005301.ont.hg38.sorted.bam SRR8955953.pacbio.hg38.sorted.bam --regions duplication_select.region.bed --out-prefix dup/duplication.pdf --fa /Files/chenjin/WGS_hg38.fa  --show-axis --show-coverage --padding 0  --show-supp
bamsnap-lrs dna --bam SRR16005301.ont.hg38.sorted.bam SRR8955953.pacbio.hg38.sorted.bam --regions duplication_select.region.bed --out-prefix dup/duplication.svg --fa /Files/chenjin/WGS_hg38.fa  --show-axis --show-coverage --padding 0  --show-supp

# RNA splicing 
bamsnap-lrs rna --bam human.hg38.Pacbio.merged.bam --regions RNA_select.region.bed --out-prefix RNA/RNA.png  --mapq 0 --fa /Files/chenjin/WGS_hg38.fa  --padding 0  --bed LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed --show-axis
bamsnap-lrs rna --bam human.hg38.Pacbio.merged.bam --regions RNA_select.region.bed --out-prefix RNA/RNA.svg  --mapq 0 --fa /Files/chenjin/WGS_hg38.fa  --padding 0  --bed LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed --show-axis
bamsnap-lrs rna --bam human.hg38.Pacbio.merged.bam --regions RNA_select.region.bed --out-prefix RNA/RNA.pdf  --mapq 0 --fa /Files/chenjin/WGS_hg38.fa  --padding 0  --bed LRGASP_manual_annotation.human.cDNA_PacBio.bamsnaplrs.bed --show-axis

# SNV phasing
bamsnap-lrs highlight --highlight-vcf HG002.phased_het_snps.vcf.gz --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam --regions phasing_select.region.bed  --out-prefix phasing/SNV.png --fa /Files/chenjin/WGS_hg38.fa  --no-hap-filter --padding 0 --mapq 0 
bamsnap-lrs highlight --highlight-vcf HG002.phased_het_snps.vcf.gz --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam --regions phasing_select.region.bed  --out-prefix phasing/SNV.svg --fa /Files/chenjin/WGS_hg38.fa  --no-hap-filter --padding 0 --mapq 0 
bamsnap-lrs highlight --highlight-vcf HG002.phased_het_snps.vcf.gz --bam HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam --regions phasing_select.region.bed  --out-prefix phasing/SNV.pdf --fa /Files/chenjin/WGS_hg38.fa  --no-hap-filter --padding 0 --mapq 0 
