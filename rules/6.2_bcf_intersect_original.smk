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
        intersection="%s/6_GATK_variants/intersection_original/0000.vcf.gz" % config["results_folder"]
    params:
        intersection="%s/6_GATK_variants/intersection_original" % config["results_folder"],
        n_files=len(lst_original)
    log:
        "%s/logs/6-2_bcftools_intersect.log" % config["results_folder"]
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
        
        bcftools isec -c none {output.filtered} -p {params.intersection} -O z --nfiles={params.n_files} &>{log}
    """
