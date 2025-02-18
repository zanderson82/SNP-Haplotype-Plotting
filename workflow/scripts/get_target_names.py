import os
import subprocess

DIR = snakemake.input.sample_dir
target_file = snakemake.input.target_file
def get_target_names(target_file):
    sample_list = []
#opens the target file and strips the lines and makes a list e.g.:
#HG00331
#HG01615
#HG01812
    with open(target_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            samples = line.strip()
            sample_list.append(samples) 
    suffix = ".PMDV_FINAL.phased.vcf.gz"
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
    
    #sample_prefix_list = snakemake.output.sample_prefix_list
    #with open(sample_prefix_list, 'w') as f:
        #for sample in sample_prefixes:
            #f.write(sample + '\n')
    
    input_vcf_list = snakemake.output.input_vcf_list
    with open(input_vcf_list, 'w') as v:
        for sample_file in new_sample_files:
            v.write(sample_file + '\n')
    
    return new_sample_files 

get_target_names(target_file)