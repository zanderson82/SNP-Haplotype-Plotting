# SNP-Haplotype-Plotting
This tool has been designed to visualize SNPs in samples with haplotype resolution.

All dependencies are included in workflow/rules/envs/

They are as follows:

* bcftools-1.19

* bedtools-2.31.1

* iPython-8.15.0

* Rtools-1.1
  - r-base=4.2.3
  - r-tidyverse
  - r-ggplot2
  - r-reshape2
  - r-readr
  - r-plyr
  - r-dplyr
  - r-optparse
  - bioconductor-biobase
  - bioconductor-biocgenerics
  - bioconductor-dnacopy

Below is a DAG that shows the flow of data through the snakemake:

![image](https://github.com/user-attachments/assets/f15d4610-166f-4bff-9d6b-ebcb703c9a2d)

The initial input is a vcf with the following structure:

![image](https://github.com/user-attachments/assets/6965693f-a640-4a07-9313-1de05abe4557)

The config.yaml file will be where you input most of the information that will change from run to run:

![image](https://github.com/user-attachments/assets/ae8e2a38-8c84-4bb2-a5f0-873b7f011349)

### The ``targetfile`` is where you will provide the sample IDs

EX: HG00151

### ``gene_target``:



### The ``sampledirectory`` will be wherever your vcfs are located

In my case they are located in /waldo/1KGP_alignments/align-card-2.24-hg38/FIRST_100

### The ``annotationfile`` and ``annotation_no_exon`` file

These two files are tab separated bed and gtf files that are already in config. As the names suggest, one of them has individual exon information while the other has just genes. To change which one you use, you will change the name in one of the rules ( Example coming up later).

### The ``outputpath`` is where the snake exists

You will change this to wherever you copy the repository to. 

### The ``demographicFile`` is a file with sample ID information and sex information. 

You can provide a file with two columns: Sample name and Sex. Down the road you can then use this to provide that information in the output plots if you are visualizing SNPs on the X chromosome. 

### The ``plot_format`` can be ``pdf``, ``svg`` or ``png``

###  The ``sample_suffix`` will change depending on how you name your files.

There is a python script that uses the sample IDs provided in the targets file and generates the ${sample} wildcards.

When taking ``HG00151`` in as the target and ``/waldo/1KGP_alignments/align-card-2.24-hg38/FIRST_100`` as the sampledirectory, it will output: ``HG00151-ONT-hg38-R9-LSK110-guppy-sup-5mC``

This output is then saved as the sample wildcard and used to find the input vcf:

${sampledirectory}/${sample}/${sample}${sample_suffix}

``/waldo/1KGP_alignments/align-card-2.24-hg38/FIRST_100/HG00151-ONT-hg38-R9-LSK110-guppy-sup-5mC/HG00151-ONT-hg38-R9-LSK110-guppy-sup-5mC.PMDV_FINAL.phased.vcf.gz``

### If you follow the ``${sampledirectory}/${sample}/${sample}${sample_suffix}`` pattern then there should be no problems finding the input vcfs

## The ``description`` will be used when making the output plot name. 
If the ``description`` is ``OPSIN-option-test-2-HG00151-SNP-PLOT``, then the output plot will be named ``OPSIN-option-test-2-HG00151-SNP-PLOT.pdf``.

These plots will be output to ``outpath/SNP-Plots/``

## Usage examples

Before running the snakemake you must make a conda environment for snakemake itself. To do this you will run:

 ``conda create -n snakemake -c conda-forge -c bioconda snakemake=8.16.0``

You must have this activated when running the tool. 

To run the tool you will run the following:

``snakemake --use-conda --conda-frontend conda --cores 1 --rerun-incomplete --rerun-triggers mtime --keep-going -p``

* ``--use-conda`` tells the snakemake to use the environments that we have told each process to use.
* ``--conda-frontend conda`` tells snakemake which package manager to use ( this is either conda or mamba )
* ``--cores`` tells snakemake how many cores to use in total. For this tool 1 is more than enough since we are not doing any resource intensive processing.
* ``--rerun-incomplete`` acts as a resume function if the snakemake runs into an error.
* ``--rerun-triggers mtime`` defines what triggers the rerunning of a job. ``mtime`` uses file modification times to determine this.
### To test you will use the command with a -n flag
This will be a dry run. It should list all of the rules that will be ran. If this works, then you will remove the n flag and proceed with the command above.

## When running for the first time, snakemake will have to install the environments provided which may take a couple minutes.

## Customization
When it comes to plotting, there are a handful of flags that can be helpful.

### Below is the ``make-SNP-plot`` rule:
```
rule make_SNP_plot:
    input:
        filtered_files = f"{config['outputpath']}/filtered-annotated-VCFs/file_list.txt"
    output:
        plots = f"{config['outputpath']}/SNP-Plots/{config['description']}-SNP-PLOT.{config['plot_format']}"
    params:
        demographic_data = config["demographicFile"],
        plotscript = "workflow/scripts/SNP-plot-with-options.R",
        title = config["description"]
    conda:
        config["Rtools-1.1"]
    shell:
        """
        Rscript {params.plotscript} -d {params.demographic_data} -o {output.plots} -f {input.filtered_files} -t {params.title} -c ""
        """
```

What matters here is the ``shell`` portion. This is where you will change the option inputs. 

### Below are the possible options that you can add or remove
```
R
option_list <- list(
  make_option(c("-d", "--demographic"), type="character", default=NULL,
              help="Path to demographic file [default= %default]"),
  make_option(c("-o", "--output"), type="character", default="plot.pdf",
              help="Output plot file [default= %default]"),
  make_option(c("-f", "--file_list"), type="character",
              help="Path to file list"),
  make_option(c("-t", "--title"), type="character", default="SNP Plot",
              help="Plot title [default= %default]"),
  make_option(c("-x", "--xchrom"), action="store_true", default=FALSE,
              help="Process as X chromosome data [default= %default]"),
  make_option(c("-c", "--gene_exons"), type="character", default=NULL,
              help="Comma-separated list of gene_exons to include (optional, empty string enables exon processing without filtering)"),  # Added comma here
  make_option(c("-w", "--width"), type="numeric", default=12,
              help="Plot width in inches [default= %default]"),
  make_option(c("-h", "--height"), type="numeric", default=8,
              help="Plot height in inches [default= %default]")
)
```

### Option usage and descriptions:
* ``-d --demographic``
  * This is where you would fill in ``{params.demographic_data}`` if you wanted to use demographic data.
  * If not, you can remove the -d flag: ``Rscript {params.plotscript} -o {output.plots} -f {input.filtered_files} -t {params.title} -c "" ``
* ``-o --output``
  * This option will always stay the same because it is connected to the ``output`` portion of the rule.
* ``-f --file_list``
  *This option will awlways be the same because the file_list is generated in a previous rule with the names of the filtered and annotated vcfs.
* ``-t --title``
  * This option will always be the same because it is referencing the ``description`` input in the config file.
* ``-x --xchrom``
  * This option, when used, will not work without a demographic file that contains the sex information.
  * When used, it will collapse the haplotype of XY individuals to one row. Only use this when you are using demographic information.
* ``-c --gene_exons``
  * This option will be passed when you choose to annotate with the gene and exon file.
  * It can be passed with a list of gene_exon names (i.e., ``-c "OPN1LW_3,OPN1LW_5"`` for the 3rd and 5th exon of the OPN1LW gene) When provided a list of names, the plot will only include variants in these exons.
  * If you have annotated with genes and exons, but want to show all that are present then you will pass ``-c ""``
* ``-w --width``
  * The width option default value is 12, however if you create the plot and decide that you want it larger, you can add the option and input your own value.
* ``-h --height``
  * The width option default value is 8, however if you create the plot and decide that you want it larger, you can add the option and input your own value.
 
#### If you want to remake a plot with the same name, you will have to delete the plot from SNP-Plots or move it somewhere else.
Snakemake works backwards to determine which steps must be ran. So if the final output is present, then it will not think to run any step. 
