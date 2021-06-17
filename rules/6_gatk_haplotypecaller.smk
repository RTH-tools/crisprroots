rule GATK_haplotypecaller:
    """
    Calling germline variants by comparing to reference
    """
    input:
        "%s/4_GATK_dedupSplit/{sample}/{sample}.Split.bam" % config["results_folder"]
    output:
        vcf="%s/6_GATK_variants/{sample}/{sample}.vcf.gz" % config["results_folder"]
    log:
        "%s/logs/6_GATK_haplotypecaller_{sample}.log" % config["results_folder"]
    params:
        reference=config["picard_reference"],
        dbsnp=config["common_variants_vcf"]
    conda:
        "../envs/gatk-picard.yaml"
    shell: """
        
        #******PARAMETERS*****
        # -R : path to genome reference
        # -I : path to bam input
        # -O : path to vcf output
        # --standard-min-confidence-threshold-for-calling : standard min confidence threshold for calling
        # --dbsnp : path to SNP annotations from dbsnp
        
        gatk HaplotypeCaller  \
        -R {params.reference} \
        -I {input} \
        -O {output.vcf} \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp {params.dbsnp} \
        &>{log}

    """
