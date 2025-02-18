import os
import subprocess

DIR = config["sampledirectory"]
target_file = config["targetfile"]
def get_target_names(target_file):
    sample_list = []
    with open(target_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            samples = line.strip()
            sample_list.append(samples) 
    suffix = config["sample_suffix"]
    new_sample_files = []
    sample_prefixes = []
    for sample in sample_list:
        find_dir_command = f"find {DIR} -maxdepth 1 -type d -name {sample}*" ## find is a shell command that looks in the specified {DIR} with a depth of 1 meaning that it will only search within the specified directory and not in any subdirectories. 
        dir_result = subprocess.getoutput(find_dir_command)
        sample_dirs = dir_result.splitlines()       
        for sample_dir in sample_dirs:
            file_prefix = sample_dir.split('/')[-1]
            new_sample_file = os.path.join(sample_dir, file_prefix + suffix)
            new_sample_files.append(new_sample_file)
            sample_prefixes.append(file_prefix)
    return sample_prefixes

rule bcftools_view:
    input:
        vcf = lambda wildcards: f"{config['sampledirectory']}/{wildcards.sample}/{wildcards.sample}{config["sample_suffix"]}"
    output:
        parsed_vcf = (f"{config['outputpath']}/parsed-VCFs/{{sample}}.parsed.vcf.gz")
    params:
        gene_coordinates = config["gene_target"]
    conda:
        "bcftools-1.19"
    shell:
        """
            bcftools view -R {params.gene_coordinates} {input.vcf} > {output.parsed_vcf}.tmp
            bgzip {output.parsed_vcf}.tmp 
            mv {output.parsed_vcf}.tmp.gz {output.parsed_vcf}
            bcftools index {output.parsed_vcf}
        """


