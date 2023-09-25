rule ROOTS_variant_screening:
    """
    Perform variant-based off-target detection and evaluation
    """
    input:
        bed="%s/10_bedops_vcf2bed/mutect2_lifted.bed" % config["results_folder"],
        dic="%s/0_utils/dict_chroms_lengths.pkl" % config["results_folder"],
        twobit="%s/0_utils/variated_genome.2bit" % config["results_folder"],
        lifted_edits="%s/0_utils/lifted_edits.bed" % config["results_folder"],
        rnafold="%s/0_utils/guide_RNA.fold" % config["results_folder"]
    output:
        tsv=temp("%s/11_VariantBasedScreening/EvaluatedVariantsOffTargets.tsv" % config["results_folder"]),
        bed=temp("%s/11_VariantBasedScreening/EvaluatedVariantsOffTargets.bed" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
        cut_position=config["Endonuclease"]["cut_position"],
        gRNA_sequence=config["Endonuclease"]["gRNA_sequence"],
        binding_sites=config["Endonuclease"]["binding_sites_seq"],
        seed_size=config["Endonuclease"]["seed_region"],
        binding_sites_ratios=config["Endonuclease"]["binding_sites_ratios"],
        binding_sites_distance=config["Endonuclease"]["binding_sites_distance"],
        extend_binding=config["Endonuclease"]["extend_binding"],
        out_folder="%s/11_VariantBasedScreening" % config["results_folder"],
        expand_search=config["VariantBasedScreening"]["expand_search"]
    log:
        "%s/logs/11_VariantScreening.log" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
        #******PARAMETERS*****
        # -v : table containing identified variants
        # -cp : number of nucleotides between the end of the guide RNA molecule and the cleavage position
        # -g : gRNA sequence
        # -b : sequence of the endonuclease binding site(s)
        # -bd : number of nucleotides between the end of the gRNA and the beginning of the endonuclease binding site
        # -t : genome in twobit format
        # -s : size of the seed region
        # --extend_binding : try to find a better binding site by extending the search on the DNA at the PAM-distal region
        # -o : path to the folder to be used for output
        # --expand_search : expand the search region for binding sites +- N nucleotides, where N is the given parameter. 
        # -p : position(s) where on-target edits are expected (bed file)
        # -d : path to dictionary of chromosome length
        # -rfo : path to output of RNAfold
        # -br : list of weights for binding site(s)
    
        mkdir -p {params.out_folder}

        python3 {params.scripts_folder}/scripts/11.0_scan_evaluate_variants_offtargets.py \
        --variants_table {input.bed} \
        --cut_position {params.cut_position} \
        --gRNA_sequence {params.gRNA_sequence} \
        --binding_sites {params.binding_sites} \
        --seed_size {params.seed_size} \
        --binding_sites_distance {params.binding_sites_distance} \
        --two_bit_genome {input.twobit} \
        --extend_binding {params.extend_binding} \
        --rna_fold_output {input.rnafold} \
        --output {params.out_folder} \
        --expand_search {params.expand_search} \
        --dict_genome_info {input.dic} \
        --binding_sites_ratios {params.binding_sites_ratios} \
        --on_target_pos {input.lifted_edits} \
        &>{log}
    """
