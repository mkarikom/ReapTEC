#!/bin/bash

### Creating a mask file used for enhancer analysis. Mask transcripts which are not lncRNAs.
### This script uses human GENCODE v39 GTF file. The same can be applied to any other GTF file from the GENCODE database. 
cd $REF_PARENT_DIR

### Remove header from GTF file
awk 'NR > 5 { print }' gencode.v41.primary_assembly.annotation.gtf > Noheader_gencode.v41.primary_assembly.annotation.gtf

### Select all transcripts from gtf file 
awk '$1 ~ /chr/ { print $0 }' Noheader_gencode.v41.primary_assembly.annotation.gtf \
| awk 'BEGIN{OFS="\t"} {if($3 == "transcript") {print $0}}' \
| sort -k1,1 -k4,4n  > Noheader_gencode.v41.chr.transcript.gtf

wc Noheader_gencode.v41.chr.transcript.gtf
### n = 244,939

### Select all non-lncRNA gene-type and transcript-type from gtf file 
awk 'BEGIN{OFS="\t"} {if($14 != "\"lncRNA\";" && $18 != "\"lncRNA\";") {print $0}}' Noheader_gencode.v41.chr.transcript.gtf > Noheader_gencode.v41.chr.gene.transcript.nonlncRNA.gtf

wc Noheader_gencode.v41.chr.gene.transcript.nonlncRNA.gtf
### n = 193,007

### GTF to bed file
awk 'BEGIN{OFS="\t"} {{print $1, $4-1, $5, $1":"$4-1".."$5","$7, "0", $7, $14"-"$18"_"$16} }' Noheader_gencode.v41.chr.gene.transcript.nonlncRNA.gtf > Noheader_gencode.v41.chr.gene.transcript.nonlncRNA.bed

### Extend 5'-end by +/-300bp
awk 'BEGIN{OFS="\t"} {if($6 == "+") {print $1, $2-300, $2+300, $1":"$2-300".."$2+300","$6, $5, $6, $7} else { print $1, $3-300, $3+300, $1":"$3-300".."$3+300","$6, $5, $6, $7} }' Noheader_gencode.v41.chr.gene.transcript.nonlncRNA.bed > Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.bed

### Sort and merge file
sort -k1,1 -k2,2n Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.bed | bedtools merge -s -c 6,7 -o distinct,collapse -delim "|" -i stdin > Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.merged.bed
wc Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.merged.bed
### n = 90,565

awk 'BEGIN{OFS="\t"} {{print $1, $2, $3, $1":"$2".."$3","$4"|"$5, "0", $4} }' Noheader_gencode.v41.chr.gene.transcript.nonlncRNA600.merged.bed > $MASK_REF