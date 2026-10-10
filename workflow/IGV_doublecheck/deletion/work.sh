# At some loci, where BamSnap-LRS did not clearly display the supporting SV reads alignment, similar limitations were observed in SVhawkeye, Wally
# DEL3	chrX	129,457,342	129,458,486
# DEL4	chr18	2,117,131	2,118,423
# DEL11	chr18	821,225	822,284
# DEL14	chr1	9,730,845	9,731,970
# DEL17	chr10	104,106,212	104,107,293
# DEL18	chrX	126,985,536	126,986,632
# DEL21	chr6	113,456,600	113,457,651
# DEL30	chr19	16,063,693	16,064,753
# DEL38	chr3	187,569,363	187,571,130
# DEL46	chr16	91,718,209	91,719,286
# Therefore, we extracted the reads from these loci and took screenshots in IGV for verification. 
# For instructions on the BAM files, please refer to https://github.com/wanyi0309/Bamsnap-LRS_visualization_gallery/blob/main/workflow/gallay_used/insertion_work.sh

samtools view -h HIFI.primary.sort.bam chrX:129457342-129458486 chr18:2117131-2118423 chr18:821225-822284 chr1:9730845-9731970 chr10:104106212-104107293 chrX:126985536-126986632 chr6:113456600-113457651 chr19:16063693-16064753 chr3:187569363-187571130 chr16:91718209-91719286 | samtools sort -o DEL.HIFI.select.bam -
samtools index DEL.HIFI.select.bam

samtools view -h ONT.primary.sort.bam chrX:129457342-129458486 chr18:2117131-2118423 chr18:821225-822284 chr1:9730845-9731970 chr10:104106212-104107293 chrX:126985536-126986632 chr6:113456600-113457651 chr19:16063693-16064753 chr3:187569363-187571130 chr16:91718209-91719286 | samtools sort -o DEL.ONT.select.bam -
samtools index DEL.ONT.select.bam
