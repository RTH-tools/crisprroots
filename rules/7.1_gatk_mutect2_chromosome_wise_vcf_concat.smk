rule GATK_mutect2_chromosome_wise_vcf_concat:
    """
    Concatenate the MuTect2 variants discovered chromosome-wise parallelly
    """
    input:
        checkpoint="%s/7_GATK_mutect2_chromosome_wise/mutect2.done" % config["results_folder"],
        Rlog="%s/0_utils/0-0_R_install_check.tmp" % config["results_folder"]
    output:
        vcfgz="%s/7_GATK_mutect2_merged/mutect2.vcf.gz" % config["results_folder"],
    log:
        "%s/logs/7_run_gatk_mutect2_chromosome_wise_vcf_concat.log" % config["results_folder"],
    params:
        scripts_folder=config["path_to_snakemake"],
        dir=directory("%s/7_GATK_mutect2_chromosome_wise/" % config["results_folder"]),
        vcf="%s/7_GATK_mutect2_merged/mutect2_merged.vcf" % config["results_folder"],
        vcf_filtered="%s/7_GATK_mutect2_merged/mutect2_filtered.vcf" % config["results_folder"],
    conda:
        "../envs/py3.yaml"
    shell: """
        
        Rscript {params.scripts_folder}/scripts/7_concatenate_vcf.R {params.dir} {params.vcf}
        cat {params.vcf} | awk '(/^#/) || ($7=="PASS")' > {params.vcf_filtered}
        bgzip -c {params.vcf_filtered} > {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        rm {params.vcf} {params.vcf_filtered}
   """
