rule BCFTOOLS_intersect_original:
    """
    Intersect germline variant calls to get variant consensus
    """
    input:
        filtered=expand("%s/6_GATK_variants/{sample}/variants_filtered.vcf" % config[
            "results_folder"],sample=lst_original)
    output:
        filtered=expand("%s/6_GATK_variants/{sample}/variants_filtered.vcf.gz" % config[
            "results_folder"],sample=lst_original),
        intersection_filtered="%s/6_GATK_variants/intersection_original/0000_filtered.vcf.gz" % config["results_folder"]
    params:
        intersection="%s/6_GATK_variants/intersection_original" % config["results_folder"],
        n_files=len(lst_original),
        scripts_folder=config["CRISPRroots"],
        intersection_filtered_vcf="%s/6_GATK_variants/intersection_original/0000_filtered.vcf" % config[
            "results_folder"]
    log:
        bcf="%s/logs/6-2_bcftools_intersect.log" % config["results_folder"],
        filter="%s/logs/6-2_filter.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
        #******PARAMETERS*****
        # -O : output type, z=compressed VCF
        # -o : path to output file

        for x in {input.filtered}; do bcftools view $x -O z -o $x.gz; done
        for x in {output.filtered}; do bcftools index $x; done

        #******PARAMETERS*****
        # -c : if none, only records with identical REF and ALT are used
        # -p : prefix dir
        # -O : output type, z=compressed VCF
        # --nfiles : --nfiles=INT the program outputs SNP positions present in INT isolates

        bcftools isec -c none {output.filtered} -p {params.intersection} -O z --nfiles={params.n_files} &>{log.bcf}
        python {params.scripts_folder}/scripts/6.2_remove_undefined_deletions.py -i {params.intersection}/0000.vcf.gz -o {params.intersection_filtered_vcf} &>{log.filter}
        bcftools view {params.intersection_filtered_vcf} -Oz -o {output.intersection_filtered}
        bcftools index {output.intersection_filtered}
    """
