rule MULTIQC_after_cleaning:
    """
    Collect reports of assessment of sequencing quality (multiqc)
    """
    input:
        "%s/preproc/1-1_fastqc_after_cutadapt_cleaning/input.lst" % config["results_folder"]
    output:
        summary_file="%s/preproc/1-2_multiqc_after_cutadapt_cleaning/multiqc_data/multiqc_general_stats.txt" % config[
            "results_folder"]
    log:
        "%s/logs/preproc/1-2_multiqc_after_cutadapt_cleaning.log" % config["results_folder"]
    params:
        folder="%s/preproc/1-2_multiqc_after_cutadapt_cleaning" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
        printf \"Producing multiQC reports after cleaning\\n\"

        #******PARAMETERS*****
        # -l : file with input files list
        # -o : path to output directory
        # -f : force re-create folder

        multiqc -l {input} -f -o {params.folder} &> {log}
    """
