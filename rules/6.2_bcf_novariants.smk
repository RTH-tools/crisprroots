rule BCFTOOLS_novariants:
    """
    Creates the output of BCFTOOLS_intersect_original (6.2_bcf_intersect_original) if variant calling is not requested
    """
    output:
        intersection="%s/6_GATK_variants/intersection_original/0000_filtered.vcf.gz" % config["results_folder"],
        vcf_intersection=temp("%s/6_GATK_variants/intersection_original/0000.vcf" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"]
    conda:
        "../envs/py3.yaml"
    shell: """
        #******PARAMETERS*****
        # -v : empty variant file        
    
        python3 {params.scripts_folder}/scripts/6.2_make_novariants_files.py \
        -v {output.vcf_intersection} 
        
        #******PARAMETERS*****
        # -O : output type, z=compressed VCF
        # -o : path to output file
        
        bcftools view {output.vcf_intersection} -Oz -o {output.intersection}
        bcftools index {output.intersection}
    """
