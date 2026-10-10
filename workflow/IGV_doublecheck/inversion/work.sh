# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# INV1	chr1	146,501,742	146,528,918
# INV17	chr5	17,508,287	17,584,034
# INV19	chr5	70,534,702	70,604,143
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/inversion_work.sh

samtools view -h HIFI.primary.sort.bam chr1:146501742-146528918 chr5:17508287-17584034 chr5:70534702-70604143 | samtools sort -o INV.HIFI.select.bam -
samtools index INV.HIFI.select.bam

samtools view -h ONT.primary.sort.bam chr1:146501742-146528918 chr5:17508287-17584034 chr5:70534702-70604143 | samtools sort -o INV.ONT.select.bam -
samtools index INV.ONT.select.bam
