rule BEDOPS_vcf2bed:
    """
    Save vcf file as a bedops file, transform it to a bed file.
    """
    input:
        vcf="%s/7-1_GATK_mutect2_merged/mutect2.vcf.gz" % config["results_folder"],
    output:
        bed=temp("%s/10_bedops_vcf2bed/mutect2_variants.bed" % config["results_folder"]),
    params:
        scripts_folder=config["CRISPRroots"],
        tmp="%s/10_bedops_vcf2bed/mutect2_full.bed" % config["results_folder"],
    log:
        bedops="%s/logs/10_bedops.log" % config["results_folder"],
    singularity: config["Singularity"]
    shell: """
        zcat {input.vcf} | vcf2bed --do-not-split --keep-header - > {params.tmp}
        python3 {params.scripts_folder}/scripts/Bedopsbed2bed.py -i {params.tmp} -o {output.bed} &>{log.bedops}
        """
