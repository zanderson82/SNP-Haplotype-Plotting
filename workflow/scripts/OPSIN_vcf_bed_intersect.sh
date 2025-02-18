#!/bin/bash

dir="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card"
demographic_file="/n/zanderson/OPSIN_BAMS/1000G_Sample_Demographic_Info.txt"
Annotated_file="/n/zanderson/3_ENSEMBL_CANONICAL_EXON_in_PROTEIN_CODING_GENE_minus_UTR_gencode.v45.annotation.gtf"
#for filtered_vcfs in $dir/*/*_BOTH_OPSINS_ALL_EXONS.vcf.gz; do
    #sample_name=$(basename "$filtered_vcfs" "_BOTH_OPSINS_ALL_EXONS.vcf.gz")
    #output_bed="${dir}/${sample_name}/${sample_name}_BOTH_OPSINS_ALL_EXONS_annotated.bed"
    #bedtools intersect -a "$filtered_vcfs" -b "$Annotated_file" -wa -wb | awk ' BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $15, $16}' > "$output_bed"
#done
for annotated_beds in $dir/*/*_BOTH_OPSINS_ALL_EXONS_annotated.bed; do
    sample_name=$(basename "$annotated_beds" "_BOTH_OPSINS_ALL_EXONS_annotated.bed")
    output_filtered_bed="${dir}/${sample_name}/${sample_name}_filtered_BOTH_OPSINS_ALL_EXONS_annotated_all_PASS_GE_5.bed"
    awk '($6 >= 5 && $7 == "PASS") || ($7 == "refCall")' $annotated_beds > $output_filtered_bed
done



