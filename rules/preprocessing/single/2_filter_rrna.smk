rule BBDUK_filter_rrna:
    """
    Cleaning of reads for ribosomal RNAs and other overrepresented RNAs (Bbduk)
    """
    input:
        "%s/preproc/1_cutadapt_cleaning/{sample}/{sample}%s" % (config["results_folder"],config["sample_suffix"]),
    output:
        temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}%s" % (config["results_folder"],config["sample_suffix"])),
    params:
        mcf=config["BBDuck"]["mcf"],
        K=config["BBDuck"]["K"],
        MAX_MEM=config["BBDuck"]["MAX_MEM"],
        ssu_rrna_silva_file=config["ssu_rrna_silva"],
        lsu_rrna_silva_file=config["lsu_rrna_silva"],
        first_ssu=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_filtered_SSU_only.fastq.gz" % config["results_folder"]),
        ssu_matched=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_matched_SSU.fastq.gz" % config["results_folder"]),
        lsu_matched=temp("%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}_matched_LSU.fastq.gz" % config["results_folder"])
    log:
        ssu="%s/logs/preproc/2_filter_rrna_SSU_{sample}.log" % config["results_folder"],
        lsu="%s/logs/preproc/2_filter_rrna_LSU_{sample}.log" % config["results_folder"]
    threads: 2
    conda:
        "../../../envs/preproc-qc.yaml"
    shell:"""

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
	bbduk.sh K={params.K} mcf={params.mcf} in={input} out={params.first_ssu} outm={params.ssu_matched} ref={params.ssu_rrna_silva_file} overwrite=t {params.MAX_MEM} &> {log.ssu}

	#remove LSU rRNA reads
	bbduk.sh K={params.K} mcf={params.mcf} in={params.first_ssu} out={output} outm={params.lsu_matched} ref={params.lsu_rrna_silva_file} overwrite=t {params.MAX_MEM} &> {log.lsu}

    """ 
