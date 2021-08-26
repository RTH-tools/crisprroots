rule CUTADAPT_cleaning:
    """
    Cleaning reads for adapter and quality filtering (cutadapt)
    """
    input:
        ["%s/{sample}%s" % (config["samples_folder"], config["sample_suffix_R1"]),
         "%s/{sample}%s" % (config["samples_folder"], config["sample_suffix_R2"])]
    output:
        first=temp("%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R1"])),
        second=temp("%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R2"]))
    params:
        phread_score=config["Cutadapt"]["phread_score"],
        adapter_file_R1=config["Cutadapt"]["adapter_R1"],
        adapter_file_R2=config["Cutadapt"]["adapter_R2"],
        pair_filter=config["Cutadapt"]["pair_filter"],
        min_length=config["Cutadapt"]["min_length"],
        other=config["Cutadapt"]["other"]
    log:
        "%s/logs/preproc/1_cutadapt.{sample}.log" % config["results_folder"]
    threads: 1
    conda:
        "../../../envs/preproc-qc.yaml"
    shell: """
        output_dir=$(dirname {output.first})
        [ ! -d \"$output_dir\" ] && mkdir -p $output_dir

        #******PARAMETERS*****
        # -q : threshold used for quality trimming
        # -a : path to file containing adapter sequence that might be ligated 3' end of the first read
        # -A : path to file containing adapter sequence that might be ligated 3' end of the second read
        # -o : path to output file for first read
        # -p : path to output file for second read
        # --minimum-length : reads shorter than this length are discarded
        # --pair-filter : if "any" then the pair is discarded if one of the reads meet the filtering criterium (min length)

        printf \"Removing low quality reads and Illumina adapter from the reads\\n\"
        cutadapt -q {params.phread_score} -a file:{params.adapter_file_R1} -A file:{params.adapter_file_R2} -o {output.first} -p {output.second} --minimum-length={params.min_length} --pair-filter={params.pair_filter} {params.other} {input} &> {log}
    """
