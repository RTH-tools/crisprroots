rule FASTQC_quality_check:
    """
    Assessment of sequencing quality (FASTQC)
    """
    input:
        ["%s/{sample}%s" % (config["samples_folder"], config["sample_suffix_R1"]),
         "%s/{sample}%s" % (config["samples_folder"], config["sample_suffix_R2"])]
    output:
        directory("%s/preproc/0_fastqc/{sample}" % config["results_folder"])
    log:
        "%s/logs/preproc/0_fastqc_{sample}.log" % config["results_folder"]
    threads: 4
    conda:
        "../../../envs/preproc-qc.yaml"
    shell: """
        [ ! -d \"{output}\" ] && mkdir {output}

        #******PARAMETERS*****
        # -t : number of threads to be used
        # -o: path to output directory

	    fastqc -t {threads} {input} -o {output} &> {log}
    """
