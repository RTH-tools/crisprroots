rule CUTADAPT_cleaning:
    """
    Cleaning reads for adapter and quality filtering (cutadapt)
    """
    input:
        "%s/{sample}%s" % (config["samples_folder"],config["sample_suffix"])
    output:
        temp("%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (config["results_folder"],config["sample_suffix"]))
    params:
        phread_score=config["Cutadapt"]["phread_score"],
        adapter_file=config["Cutadapt"]["adapter"],
        min_length=config["Cutadapt"]["min_length"]
    log:
        "%s/logs/preproc/1_cutadapt_{sample}_log" % config["results_folder"]
    threads: 1
    conda:
        "../../../envs/preproc-qc.yaml"
    shell:"""

        output_dir=$(dirname {output})
        [ ! -d \"$output_dir\" ] && mkdir -p $output_dir

        #******PARAMETERS*****
        # -q : threshold used for quality trimming
        # -a : path to file containing adapter sequence that might be ligated 3' end of the first read
        # -o : path to output file for first read
        # --minimum-length : reads shorter than this length are discarded

        printf \"Removing low quality reads and Illumina adapter from the reads\\n\"
        cutadapt -q {params.phread_score} -a file:{params.adapter_file} -o {output} --minimum-length={params.min_length} {input} &> {log}
     """
