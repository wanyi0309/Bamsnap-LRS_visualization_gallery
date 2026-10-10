# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# INS7	chr16	94,802,848	94,803,148
# INS8	chr8	56,288,249	56,288,549
# INS13	chr14	48,411,966	48,412,266
# INS16	chr7	1,039,175	1,039,475
# INS19	chr6	55,938,971	55,939,271
# INS21	chr6	169,578,532	169,578,832
# INS29	chr17	12,350,231	12,350,531
# INS31	chr1	1,558,949	1,559,249
# INS37	chr3	21,154,107	21,154,407
# INS42	chr20	65,168,891	65,169,191
# INS44	chr1	104,718,432	104,718,732
# INS45	chr3	178,976,521	178,976,821
# INS46	chr17	83,460,369	83,460,669
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/insertion_work.sh

samtools view -h HIFI.primary.sort.bam chr16:94802848-94803148 chr8:56288249-56288549 chr14:48411966-48412266 chr7:1039175-1039475 chr6:55938971-55939271 chr6:169578532-169578832 chr17:12350231-12350531 chr1:1558949-1559249 chr3:21154107-21154407 chr20:65168891-65169191 chr1:104718432-104718732 chr3:178976521-178976821 chr17:83460369-83460669 | samtools sort -o INS.HIFI.select.bam -
samtools index INS.HIFI.select.bam

samtools view -h ONT.primary.sort.bam chr16:94802848-94803148 chr8:56288249-56288549 chr14:48411966-48412266 chr7:1039175-1039475 chr6:55938971-55939271 chr6:169578532-169578832 chr17:12350231-12350531 chr1:1558949-1559249 chr3:21154107-21154407 chr20:65168891-65169191 chr1:104718432-104718732 chr3:178976521-178976821 chr17:83460369-83460669 | samtools sort -o INS.ONT.select.bam -
samtools index INS.ONT.select.bam
