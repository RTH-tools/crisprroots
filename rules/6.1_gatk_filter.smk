rule GATK_filters_haplotypecaller:
    """
    Filter germline variant calls
    """
    input:
        vcf="%s/6_GATK_variants/{sample}/{sample}.vcf.gz" % config["results_folder"]
    output:
        filters=temp("%s/6_GATK_variants/{sample}/{sample}.filters.vcf.gz" % config["results_folder"]),
        filtered="%s/6_GATK_variants/{sample}/variants_filtered.vcf" % config["results_folder"]
    log:
        "%s/logs/6-1.GATK_filter_{sample}.log" % config["results_folder"]
    params:
        reference=config["picard_reference"],
    conda:
        "../envs/gatk-picard.yaml"
    shell: """
        #Filter variants
        # QD : quality calls 
        # FS : strand bias 
        # DP : read depth
        # cluster : number of SNPs that constitute a cluster
        # window : window-size in which SNP clusters are evaluated 
        
        gatk VariantFiltration \
        -R {params.reference} \
        -V {input.vcf} \
        -window 35 -cluster 3 \
        --filter-name FS --filter-expression "FS > 30.0" \
        --filter-name QD --filter-expression "QD < 2.0" \
        --filter-name DP --filter-expression "DP < 10" \
        -O {output.filters} \
        &>{log}

        zcat {output.filters} | awk '(/^#/) || ($7=="PASS"&&$6>200)' > {output.filtered}
    """
