rule eSNPKaryotyping_launch:
    """
    Use eSNPKaryotyping to assess chromosome integrity
    """
    input:
        vcf_edited="%s/6_GATK_variants/{sample}/variants_filtered.vcf" % config["results_folder"],
        bam_edited="%s/2_sortaligned/{sample}/Aligned.Sorted.bam" % config["results_folder"],
        Rlog="%s/0_utils/0-0_R_install_check.tmp" % config["results_folder"]
    output:
        out_image="%s/eSNPKaryotyping/{sample}_Zygosity_Blocks.pdf" % config["report_folder"]
    params:
        organism=config["species"],
        picard_ref_dict=config["picard_reference"][:-3] + ".dict",
        scripts_folder=config["path_to_snakemake"],
        report=config["report_folder"],
        dbSNP_142_per_chromosome=config["dbSNP142"]
    log:
        "%s/logs/14_eSNPKaryotype_{sample}.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # (1) : Directory containing vairants
        # (2) : Scientific name of the organism
        # (3) : Path to directory with dbSNP files
        # (4) : Path to reference genome
        # (5) : Path to bam files
        # (6) : Output folder and prefix for the output file
        
        vcf_edited_dir=`dirname {input.vcf_edited}`
        bam_edited_dir=`dirname {input.bam_edited}`
        Rscript {params.scripts_folder}/scripts/14_eSNPKaryotyping.R \
        $vcf_edited_dir \
        {params.organism} \
        {params.dbSNP_142_per_chromosome} \
        {params.picard_ref_dict} \
        $bam_edited_dir \
        {params.report}/eSNPKaryotyping/{wildcards.sample} \
        &>{log}
    """
