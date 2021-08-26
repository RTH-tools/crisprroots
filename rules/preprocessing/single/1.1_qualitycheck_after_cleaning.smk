rule FASTQC_qualitycheck_after_cleaning:
    """
    Assessment of sequencing quality (FASTQC) of cleaned reads
    """
    input:
        "%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (config["results_folder"],config["sample_suffix"])
    output:
        directory("%s/preproc/1-1_fastqc_after_cutadapt_cleaning/{sample}" % config["results_folder"])
    log:
        "%s/logs/preproc/1-1_fastqc_after_cutadapt_cleaning.{sample}.log" % config["results_folder"]
    threads: 4
    conda:
        "../../../envs/preproc-qc.yaml"
    shell:"""

        [ ! -d \"{output}\" ] && mkdir {output}

        #******PARAMETERS*****
        # -t : number of threads to be used
	    # -o: path to output directory

        fastqc -t {threads} {input} -o {output} &>{log}
    """
