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
        Rscript {params.plotscript} -o {output.plots} -f {input.filtered_files} -title {params.title} 
        """