# At three loci (B3GAT2, DEGS2, and PEX12), BamSnap-LRS displayed only a subset of the isoforms
# B3GAT2	ENSG00000112309	chr6	70,853,482	70,957,054
# DEGS2	ENSG00000168350	chr14	100,143,980	100,192,849
# PEX12	ENSG00000108733	chr17	35,574,794	35,578,963
# Therefore, we extracted the reads from these three loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/RNA_work.sh

samtools view -bh human.hg38.Pacbio.unfilter.merged.bam  chr6:70,852,482-70,958,054 chr14:100,142,980-100,193,849 chr17:35,573,794-35,579,963 > problem.gene.bam
samtools index problem.gene.bam
