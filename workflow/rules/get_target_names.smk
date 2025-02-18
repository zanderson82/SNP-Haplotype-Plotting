rule get_target_names:
    input: 
        sample_dir = config["sampledirectory"],
        target_file = config["targetfile"]
    output:
        input_vcf_list = temp("/n/zanderson/SNP-plot-snake-testing-Waldo/config/input_vcf_list.txt")
    conda:
        config["iPython-8.15.0"]
    script:
        "/n/zanderson/SNP-plot-snake-testing-Waldo/workflow/scripts/get_target_names.py"

