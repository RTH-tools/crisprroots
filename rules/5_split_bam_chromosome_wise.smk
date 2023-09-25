rule ROOTS_split_bam_chromosome_wise:
    """
    Use samtools to split bam chromosome-wise
    """
    input:
        bam="%s/4_GATK_dedupSplit/{sample}/{sample}.Split.bam" % config["results_folder"],
        idx="%s/4_GATK_dedupSplit/{sample}/{sample}.Split.bai" % config["results_folder"],
    output:
        outdir=directory("%s/5_chromosome_wise_splitted_bam/{sample}/" % config["results_folder"]),
        chrsplit_checkpoint="%s/5_chromosome_wise_splitted_bam/{sample}/split.done" % config["results_folder"]
    log:
        "%s/logs/5_split_bam_chromosome_wise_unique_{sample}.log" % config["results_folder"],
    params:
        scripts_folder=config["CRISPRroots"],
        seq_type=config["sequencing"],
        dir=directory("%s/4_GATK_dedupSplit/" % config["results_folder"])
    singularity: config["Singularity"]
    shell: """
        mkdir -p {output.outdir}
        
        #******PARAMETERS*****
        # (1) : Path to out folder
        # (2) : Sequencing method (single, paired)
        # (3) : Path to input folder where bams are located
        
        # list the chr to which reads are assigned
        {params.scripts_folder}/scripts/5_chr_to_split.sh "{output.outdir}" "{params.seq_type}" {params.dir} &>{log}
        touch "{output.chrsplit_checkpoint}"

    """
