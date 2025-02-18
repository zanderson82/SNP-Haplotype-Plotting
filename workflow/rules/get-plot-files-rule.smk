import os
import subprocess

DIR = config["sampledirectory"]
target_file = config["targetfile"]
def get_sample_prefixes(target_file):
    sample_list = []
    with open(target_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            samples = line.strip()
            sample_list.append(samples) 
    sample_prefixes = []
    for sample in sample_list:
        find_dir_command = f"find {DIR} -maxdepth 1 -type d -name {sample}*"
        dir_result = subprocess.getoutput(find_dir_command)
        sample_dirs = dir_result.splitlines()       
        for sample_dir in sample_dirs:
            file_prefix = sample_dir.split('/')[-1]
            sample_prefixes.append(file_prefix)
    return sample_prefixes

rule get_plot_files_rule:
    input:
        samples = expand(f"{config['outputpath']}/filtered-annotated-VCFs/{{sample}}.annotated.filtered.bed", sample=get_sample_prefixes(config["targetfile"]))
    output:
        plotted_files = f"{config['outputpath']}/filtered-annotated-VCFs/file_list.txt"
    conda:
        "iPython-8.15.0"
    script:
        "/n/zanderson/SNP-plot-snake-testing/workflow/scripts/get-plot-files.py"