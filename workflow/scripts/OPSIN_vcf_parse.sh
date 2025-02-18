#!/bin/bash
dir="/n/1000g/align-card-2.24-hg38/FIRST_100/"
OPN1LW="/n/zanderson/OPSIN_BAMS/OPN1LW_exons_3.5_coords.txt"
OPN1MW="/n/zanderson/OPSIN_BAMS/OPN1MW_exons_3.5_coords.txt"
BOTH_OPSINS_all_exons="/n/zanderson/OPSIN_BAMS/OPSIN_exon_coords.txt"
for vcf in $dir/*/*PMDV_FINAL.phased.vcf.gz ; do
    sample_name=$(basename "$vcf" | cut -d '-' -f 1)
output_vcf="/n/zanderson/OPSIN_BAMS/1000g_OPSIN_BAMS_card/${sample_name}/${sample_name}_BOTH_OPSINS_ALL_EXONS.vcf.gz"
    echo "output vcf: $output_vcf"
    bcftools view -R "$BOTH_OPSINS_all_exons" "$vcf" | bgzip > "$output_vcf"
    if [ $? -ne 0 ]; then
        echo "bcftools vew failed for $vcf"
        continue
    fi
    bcftools index "$output_vcf"
    ls $output_vcf
done

