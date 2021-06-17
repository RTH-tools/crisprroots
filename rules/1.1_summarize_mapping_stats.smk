rule ROOTS_summarize_mapping_stats:
    """
        Creates a table that summarizes the mapping statistics
    """
    input:
        expand("%s/1_star_align2pass/{sample}/{sample}.Log.final.out" % config["results_folder"],sample=lst_samples)
    output:
        "%s/../report/mapping_stats.xlsx" % config["results_folder"]
    log:
        "%s/logs/1-1_mapping_stats.log" % config["results_folder"]
    params:
        scripts_folder=config["path_to_snakemake"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -i : Path to stats output from star
        # -o : Path to the output summary file
    
        python3 {params.scripts_folder}/scripts/1.1_gather_mapping_stats.py \
        -i {input} -o {output} &>{log}
    """
