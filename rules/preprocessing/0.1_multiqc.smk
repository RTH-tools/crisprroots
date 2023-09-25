rule MULTIQC_input:
    """
    Collect reports of assessment of sequencing quality (multiqc)
    """
    input:
        "%s/preproc/0_fastqc/input.lst" % config["results_folder"]
    output:
        "%s/preproc/0-1_multiqc/multiqc_data/multiqc_general_stats.txt" % config["results_folder"]
    params:
        folder="%s/preproc/0-1_multiqc" % config["results_folder"]
    log:
        "%s/logs/preproc/0-1_multiqc.log" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
        printf \"Producing multiQC reports\\n\"

        #******PARAMETERS*****
        # -l : file with input files list
        # -o : path to output directory
        # -f : force re-create folder
        # rm -r {params.folder}/multiqc_data

        multiqc -l {input} -f -o {params.folder} &> {log}
    """
