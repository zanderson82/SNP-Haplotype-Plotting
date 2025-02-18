rule vcf_filter:
    input:
        annotated_vcf = f"{config['outputpath']}/annotated-VCFs/{{sample}}.annotated.vcf"
    output:
        filtered_vcf = f"{config['outputpath']}/filtered-annotated-VCFs/{{sample}}.annotated.filtered.bed"
    shell:
        """
        filtered_vcf_temp="{output.filtered_vcf}.temp"
        awk 'BEGIN {{OFS="\t"}} ($6 >= 5 && $7 == "PASS") || ($7 == "refCall")' {input.annotated_vcf} > $filtered_vcf_temp
        (echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tGENE\tEXON"; cat $filtered_vcf_temp) > {output.filtered_vcf}

        rm $filtered_vcf_temp
        """