rule get_sample_prefixes:
    input:
        sample_dir = config["sampledirectory"],
        target_file = config["targetfile"]
    output:
        sample_prefix_list = "/n/zanderson/SNP-plot-snake-testing-Waldo/config/sample_prefix_list.txt"
    conda:
        "iPython-8.15.0"
    script:
        "/n/zanderson/SNP-plot-snake-testing-Waldo/workflow/scripts/get_sample_prefixes.py"