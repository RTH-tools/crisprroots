rule ROOTS_ontarget_check_knockout:
    """
    Check for the presence of transcript KO
    """
    input:
        diffexp="%s/12-1_DESeq2/diffexpr-results.tsv" % config["results_folder"]
    output:
        report="%s/../report/on_target_knockout.xlsx" % config["results_folder"]
    log:
        "%s/logs/8_ROOTS_ontarget_check.log" % config["results_folder"]
    params:
        scripts_folder=config["CRISPRroots"],
        targets=config["Edits"]["KO"],
    singularity: config["Singularity"]
    shell: """
    
        #******PARAMETERS*****
        #-d : path to table with differentially expressed genes
        #-o : path to output report about the knockout
        #-KO : list of knockout genes
    
        python3 {params.scripts_folder}/scripts/8_check_KO.py \
                -d {input.diffexp} \
                -o {output.report} \
                -KO {params.targets} \
                &>{log}
        """
