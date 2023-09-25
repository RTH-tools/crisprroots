rule MULTIQC_after_rRNA_removal:
    """
    Collect reports of assessment of sequencing quality (multiqc)
    """
    input:
        "%s/preproc/2-1_fastqc_after_rRNA_removal/input.lst" % config["results_folder"]
    output:
        "%s/preproc/2-2_multiqc_after_rRNA_removal/multiqc_data/multiqc_general_stats.txt" % config["results_folder"]
    log:
        "%s/logs/preproc/2-2_multiqc_after_rRNA_removal.log" % config["results_folder"]
    params:
        folder="%s/preproc/2-2_multiqc_after_rRNA_removal" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
        printf \"Producing multiQC reports after cleaning\\n\"

        #******PARAMETERS*****
        # -l : file with input files list
        # -o : path to output directory
        # -f : force re-create folder

        multiqc -l {input} -f -o {params.folder} &> {log}
    """
