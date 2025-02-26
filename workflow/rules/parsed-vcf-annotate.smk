rule vcf_annotate:
    input:
        parsed_vcf=f"{config['outputpath']}/parsed-VCFs/{{sample}}.parsed.vcf.gz"
    output:
        annotated_vcf=f"{config['outputpath']}/annotated-VCFs/{{sample}}.annotated.vcf"
    params:
        annotation_file=config["annotation_no_exon"],
        outputpath=config["outputpath"]
    conda:    
        config["bedtools-2.31.1"]
    shell:
        """
        unzipped_vcf="{params.outputpath}/parsed-VCFs/{wildcards.sample}.unzipped.vcf"
        zcat {input.parsed_vcf} > $unzipped_vcf
        if [ ! -f $unzipped_vcf ]; then
            echo "Error: $unzipped_vcf not created"
            exit 1
        fi
        bedtools intersect -a $unzipped_vcf -b {params.annotation_file} -wa -wb | awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14}}' > {output.annotated_vcf} 
        rm $unzipped_vcf

        ## without exons {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14}}'
        ## with exons bedtools intersect -a $unzipped_vcf -b {params.annotation_file} -wa -wb | awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $15, $16}}' > {output.annotated_vcf} 

        """