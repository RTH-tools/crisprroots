rule ROOTS_make_utils:
    input:
        genome="%s/6_GATK_variants/variated_genome.fa" % config["results_folder"],
        chain_liftover="%s/6_GATK_variants/chain_liftover_variated_genome.chain" % config["results_folder"],
        edits="%s/0_utils/edits.bed" % config["results_folder"]
    output:
        dic="%s/0_utils/dict_chroms_lengths.pkl" % config["results_folder"],
        twobitgenome=temp("%s/0_utils/variated_genome.2bit" % config["results_folder"]),
        lifted_annotations="%s/0_utils/lifted_annotations.gtf" % config["results_folder"],
        unlifted_annotations="%s/0_utils/unlifted_annotations.gtf" % config["results_folder"],
        lifted_edits="%s/0_utils/lifted_edits.bed" % config["results_folder"],
        unlifted_edits="%s/0_utils/unlifted_edits.bed" % config["results_folder"],
        tmp=temp("%s/0_utils/tmp.bed" % config["results_folder"]),
        SNPdb_bed=temp("%s/0_utils/SNPdb.bed" % config["results_folder"]),
        lifted_snpdb="%s/0_utils/lifted_SNPdb.bed" % config["results_folder"],
        unlifted_snpdb="%s/0_utils/unlifted_SNPdb.bed" % config["results_folder"],
        lifted_rpkm="%s/0_utils/lifted_repeatmask.bed" % config["results_folder"],
        unlifted_rpkm="%s/0_utils/unlifted_repeatmask.bed" % config["results_folder"]
    log:
        dict="%s/logs/0_chrom_dict.log" % config["results_folder"],
        twobit="%s/logs/0_twobit_genome.log" % config["results_folder"],
        annot="%s/logs/0_lifted_annotations.log" % config["results_folder"],
        edits="%s/logs/0_lifted_edits.log" % config["results_folder"],
        snpdb="%s/logs/0_lifted_snpdb.log" % config["results_folder"],
        rpkm="%s/logs/0_lifted_repeatmask.log" % config["results_folder"]
    conda:
        "../envs/twobit.yaml"
    params:
        scripts_folder=config["path_to_snakemake"],
        min_match=config["Liftover"]["min_match"],
        annot=config["annotations_gtf"],
        SNPdb=config["common_variants_vcf"],
        rpkm=config["repeatmasked_regions"]
    shell: """
    
        #******PARAMETERS*****
        # liftover:
            # --minMatch : Min ratio of bases that must remap
        
        faToTwoBit {input.genome} {output.twobitgenome} &>{log.twobit}

        python3 {params.scripts_folder}/scripts/Twobitinfo.py {output.twobitgenome} {output.dic} &>{log.dict}

        liftOver \
        -minMatch={params.min_match} \
        -gff \
        {params.annot} \
        {input.chain_liftover} \
        {output.lifted_annotations} \
        {output.unlifted_annotations} &>{log.annot}

        liftOver \
        -minMatch={params.min_match} \
        {input.edits} \
        {input.chain_liftover} \
        {output.lifted_edits} \
        {output.unlifted_edits} &>{log.edits}

        zcat {params.SNPdb} | vcf2bed --do-not-split >{output.tmp}
        awk -v OFS='\t' '{{print $1,$2,$3}}' {output.tmp} >{output.SNPdb_bed}

        liftOver \
        -minMatch={params.min_match} \
        {output.SNPdb_bed} \
        {input.chain_liftover} \
        {output.lifted_snpdb} \
        {output.unlifted_snpdb} &>{log.snpdb}

        liftOver \
        -minMatch={params.min_match} \
        {params.rpkm} \
        {input.chain_liftover} \
        {output.lifted_rpkm} \
        {output.unlifted_rpkm} &>{log.rpkm}
    """
