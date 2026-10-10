# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# DUP39	chr14	105,626,058	105,630,340
# DUP47	chr16	34,570,634	34,583,930
# DUP48	chr16	34,580,005	34,587,612
# DUP49	chr16	34,590,756	34,596,605
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/insertion_work.sh

samtools view -h -M SRR16005301.ont.hg38.sorted.bam chr14:105626058-105630340 chr16:34570634-34583930 chr16:34580005-34587612 chr16:34590756-34596605 | samtools sort -o DUP.hg38.ONT.select.bam -
samtools index DUP.hg38.ONT.select.bam

samtools view -h -M SRR8955953.pacbio.hg38.sorted.bam chr14:105626058-105630340 chr16:34570634-34583930 chr16:34580005-34587612 chr16:34590756-34596605 | samtools sort -o DUP.hg38.Pacbio.select.bam -
samtools index DUP.hg38.Pacbio.select.bam
