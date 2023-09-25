rule FASTQC_qualitycheck_after_removal_rRNA:
    """
    Assessment of sequencing quality (FASTQC) of cleaned reads after removing rRNA reads
    """
    input:
        "%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}%s" % (config["results_folder"],config["sample_suffix"])
    output:
        directory("%s/preproc/2-1_fastqc_after_rRNA_removal/{sample}" % config["results_folder"])
    log:
        "%s/logs/preproc/2-1_fastqc_after_rRNA_removal_{sample}.log" % config["results_folder"]
    threads: 4
    singularity: config["Singularity"]
    shell:"""
    
        [ ! -d \"{output}\" ] && mkdir {output}

        #******PARAMETERS*****
        # -t : number of threads to be used
        # -o: path to output directory

        fastqc -t {threads} {input} -o {output} &>{log}
    """
