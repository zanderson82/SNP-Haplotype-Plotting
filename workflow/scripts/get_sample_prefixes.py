import os
import subprocess

DIR = snakemake.input.sample_dir
target_file = snakemake.input.target_file
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
    
    sample_prefix_list = snakemake.output.sample_prefix_list
    with open(sample_prefix_list, 'w') as f:
        for sample in sample_prefixes:
            f.write(sample + '\n')
    return sample_prefixes