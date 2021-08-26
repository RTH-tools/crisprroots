rule ROOTS_split_bam_chromosome_wise:
    """
    Use samtools to split bam chromosome-wise
    """
    input:
        edited=expand("%s/4_GATK_dedupSplit/{edited}/{edited}.Split.bam" % config["results_folder"],edited=lst_edited),
        original=expand("%s/4_GATK_dedupSplit/{original}/{original}.Split.bam" % config[
            "results_folder"],original=lst_original),
        edited_idx=expand("%s/4_GATK_dedupSplit/{edited}/{edited}.Split.bai" % config["results_folder"],edited=lst_edited),
        original_idx=expand("%s/4_GATK_dedupSplit/{original}/{original}.Split.bai" % config[
            "results_folder"],original=lst_original)
    output:
        chrsplit_edited=directory(expand("%s/5_chromosome_wise_splitted_bam/{edited}/" % config[
            "results_folder"],edited=lst_edited)),
        chrsplit_original=directory(expand("%s/5_chromosome_wise_splitted_bam/{original}/" % config[
            "results_folder"],original=lst_original)),
        chrsplit_checkpoint=temp("%s/5_chromosome_wise_splitted_bam/split.done" % config["results_folder"])
    log:
        edited="%s/logs/5_split_bam_chromosome_wise_unique_edited.log" % config["results_folder"],
        original="%s/logs/5_split_bam_chromosome_wise_unique_original.log" % config["results_folder"],
    params:
        scripts_folder=config["CRISPRroots"],
        seq_type=config["sequencing"],
        edited_dir=directory("%s/4_GATK_dedupSplit/" % config["results_folder"]),
        original_dir=directory("%s/4_GATK_dedupSplit/" % config["results_folder"]),
    conda:
        "../envs/py3.yaml"
    shell: """
        mkdir -p {output.chrsplit_edited}
        mkdir -p {output.chrsplit_original}
        
        #******PARAMETERS*****
        # (1) : Path to out folder
        # (2) : Sequencing method (single, paired)
        # (3) : Path to input folder where bams are located
        
        # list the chr to which reads are assigned
        {params.scripts_folder}/scripts/5_chr_to_split.sh "{output.chrsplit_edited}" "{params.seq_type}" {params.edited_dir} &>{log.edited}
        {params.scripts_folder}/scripts/5_chr_to_split.sh "{output.chrsplit_original}" "{params.seq_type}" {params.original_dir} &>{log.original}
        touch "{output.chrsplit_checkpoint}"

    """
