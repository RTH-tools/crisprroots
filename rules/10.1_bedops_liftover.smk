rule BEDOPS_liftover:
    """
    Lift variants to the coordinates of the variated genome.
    """
    input:
        bed="%s/10_bedops_vcf2bed/mutect2_variants.bed" % config["results_folder"],
        chain_liftover="%s/6_GATK_variants/chain_liftover_variated_genome.chain" % config["results_folder"]
    output:
        lifted_coords=temp("%s/10_bedops_vcf2bed/mutect2_lifted.bed" % config["results_folder"]),
        unlifted_coords=temp("%s/10_bedops_vcf2bed/mutect2_unlifted.bed" % config["results_folder"])
    params:
        min_match=config["Liftover"]["min_match"]
    log:
        "%s/logs/10-1_liftover.log" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
    
        #******PARAMETERS*****
        # liftover:
            # --minMatch : Min ratio of bases that must remap
            
        liftOver \
        -minMatch={params.min_match} \
        {input.bed} \
        {input.chain_liftover} \
        {output.lifted_coords} \
        {output.unlifted_coords} &>{log} || touch {output.lifted_coords} {output.unlifted_coords}
    """
