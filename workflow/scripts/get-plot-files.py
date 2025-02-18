import sys 
samples = snakemake.input.samples
output_file = snakemake.output.plotted_files

def main(samples, output_file):
    with open(output_file, 'w') as f:
        for sample in samples:
            f.write(sample + '\n')


main(samples, output_file)
