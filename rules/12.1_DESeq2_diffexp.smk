rule DESEQ2_differential_expression:
    """
    Differential expression analysis (DESeq2)
    """
    input:
        featurecounts="%s/9_featurecounts_quantification/feature_count.tsv" % config["results_folder"],
        Rlog="%s/0_utils/0-0_R_install_check.tmp" % config["results_folder"]
    output:
        diffexp_table="%s/12-1_DESEq2DE/diffexpr-results.tsv" % config["results_folder"]
    params:
        phenotypes=config["samples_table"],
        design_formula=config["DESEq2"]["formula"],
        scripts_folder=config["path_to_snakemake"],
        KO_genes=config["Edits"]["KO"],
        diffexp="%s/12-1_DESEq2DE" % config["results_folder"],
        report=config["report_folder"]
    log:
        s="%s/logs/12-1_setup_de.log" % config["results_folder"],
        r="%s/logs/12-1_diffexpr_r.log" % config["results_folder"],
        p="%s/logs/12-1_diffexpr_summarize.log" % config["results_folder"]
    threads: 16
    conda:
        "../envs/py3.yaml"
    shell: """
        
        #******PARAMETERS*****
        # 12.1_deseq2.R : 
            # (1) : path to workdir
            # (2) : path to samples phenotypes table
            # (3) : path to featurecount output
            # (4) : design formula
        # 12.1_find_DE_genes.py:
            # -o : path to diffexp folder
            # -KO : IDs of the KO genes
            # -d : path to diffexp table
            # -p : path to samples phenotypes table
            # -r : path to report folder
       
        {params.scripts_folder}/scripts/12.1_setup_differential_expression.py -i {input.featurecounts} -o {params.diffexp} &>{log.s}
        Rscript {params.scripts_folder}/scripts/12.1_deseq2.R {params.diffexp} {params.phenotypes} {params.diffexp}/feature_count.tsv \"{params.design_formula}\" &>{log.r}
        {params.scripts_folder}/scripts/12.1_find_DE_genes.py -o {params.diffexp} -KO {params.KO_genes} -d {output.diffexp_table} -p {params.phenotypes} -r {params.report} &>{log.p}
"""
