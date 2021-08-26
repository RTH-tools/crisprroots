rule BBDUK_filter_rrna:
    """
    Cleaning of reads for ribosomal RNAs and other overrepresented RNAs (Bbduk)
    """
    input:
        first="%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R1"]),
        second="%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R2"])
    output:
        first=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R1"])),
        second=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}%s" % (
            config["results_folder"], config["sample_suffix_R2"])),
        first_ssu=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_filtered_SSU_only%s" % (
        config["results_folder"],config["sample_suffix_R1"])),
        second_ssu=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_filtered_SSU_only%s" % (
            config["results_folder"], config["sample_suffix_R2"])),
        ssu_matched=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_matched_SSU.fastq.gz" % config[
            "results_folder"]),
        lsu_matched=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_matched_LSU.fastq.gz" % config[
            "results_folder"])
    params:
        mcf=config["BBDuck"]["mcf"],
        K=config["BBDuck"]["K"],
        MAX_MEM=config["BBDuck"]["MAX_MEM"],
        ssu_rrna_silva_file=config["ssu_rrna"],
        lsu_rrna_silva_file=config["lsu_rrna"],
    log:
        ssu="%s/logs/preproc/2_filter_rrna_SSU_{sample}.log" % config["results_folder"],
        lsu="%s/logs/preproc/2_filter_rrna_LSU_{sample}.log" % config["results_folder"]
    threads: 2
    conda:
        "../../../envs/preproc-qc.yaml"
    shell: """

        #********PARAMETERS*******
        # K : size of the Kmers used by BBduk. BBduk supports K between 1 and 31. Kmer size >31 (eg.40) implies that a read matches with a reference kmer if it matches with kmers of size 31 for position 0-31,1-32,2-33...9-40
        # mcf : mcf=0.5 requires 50% of the read bases to be covered by reference kmers
        # MAX_MEM : max amount of memory in GB to be used. Default behavior is to autodetect the available memory and use half of it, but given a large reference dataset we need to use more memory (80%ca.).
        # in1, in2 : are the paired reads given in input (output reads of cutadapt)
        # out1, out2 : are the reads that did not match any of the filters
        # outm1, outm2 : are the reads that matched with a filter (if any of the fw-reverse reads fail, both will go here)
        # ref : is the reference dataset of rRNA

        printf \"Filtering rRNAs with BBDUK\\n\"

        #remove SSU rRNA reads
        bbduk.sh K={params.K} mcf={params.mcf} in1={input.first} in2={input.second} out1={output.first_ssu} out2={output.second_ssu} outm={output.ssu_matched} ref={params.ssu_rrna_silva_file} overwrite=t {params.MAX_MEM} &> {log.ssu}

        #remove LSU rRNA reads
        bbduk.sh K={params.K} mcf={params.mcf} in1={output.first_ssu} in2={output.second_ssu} out1={output.first} out2={output.second} outm={output.lsu_matched} ref={params.lsu_rrna_silva_file} overwrite=t {params.MAX_MEM} &> {log.lsu}
    """
