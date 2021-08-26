rule ROOTS_expressionbasedscreening:
    """
    Perform expression-based off-target detection and evaluation
    """
    input:
        collapsed_coords="%s/12-4_CollapseCoordinatesGenesOff/collapsed_genes_offtargets_coordinates_singlechr_filtered.bed" %
                         config["results_folder"],
        dic="%s/0_utils/dict_chroms_lengths.pkl" % config["results_folder"],
        two_bit_genome="%s/0_utils/variated_genome.2bit" % config["results_folder"],
        rnafold="%s/0_utils/guide_RNA.fold" % config["results_folder"]
    output:
        tsv=temp("%s/12-5_ExpressionBasedScreening/EvaluatedExpressionOffTargets.tsv" % config["results_folder"]),
        bed=temp("%s/12-5_ExpressionBasedScreening/EvaluatedExpressionOffTargets.bed" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
        cut_position=config["Endonuclease"]["cut_position"],
        gRNA_sequence=config["Endonuclease"]["gRNA_sequence"],
        binding_sites=config["Endonuclease"]["binding_sites_seq"],
        seed_size=config["Endonuclease"]["seed_region"],
        binding_sites_ratios=config["Endonuclease"]["binding_sites_ratios"],
        binding_sites_distance=config["Endonuclease"]["binding_sites_distance"],
        extend_binding=config["Endonuclease"]["extend_binding"],
        out_folder="%s/12-5_ExpressionBasedScreening" % config["results_folder"]
    log:
        "%s/logs/12-5_ExpressionScreening.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
        #******PARAMETERS*****
        # -input : table containing predicted off-target coordinates intersected with gene differential expression data
        # -gRNA_sequence : gRNA sequence
        # -binding_sites_distance : Number of nucleotides between the end of the gRNA and the beginning of the endonuclease binding site
        # -seed_size : Size of the seed region
        # -two_bit_genome : Genome in twobit format
        # --extend_binding : Try to find a better binding site by extending the search on the DNA at the PAM-distal region
        # -output : Path to the folder to be used for output
        # -cut_position : Number of nucleotides from the 3p end of the predicted off-target coordinates to the cleavage position
        # --crisproff : Flag. Use if the input is the crisproff tsv output table
        # dict_genome_info : Path to dictionary of chromosome length
        # rna_fold_output : Path to output of RNAfold
        # binding_sites : Sequence of the endonuclease binding site(s)
        # binding_sites_ratios : List of weights for binding site(s)

        python3 {params.scripts_folder}/scripts/12.6_scan_eval_expression_offtargets.py \
        --input {input.collapsed_coords} \
        --cut_position {params.cut_position} \
        --gRNA_sequence {params.gRNA_sequence} \
        --seed_size {params.seed_size} \
        --two_bit_genome {input.two_bit_genome} \
        --extend_binding {params.extend_binding} \
        --rna_fold_output {input.rnafold} \
        --binding_sites {params.binding_sites} \
        --binding_sites_ratios {params.binding_sites_ratios} \
        --binding_sites_distance {params.binding_sites_distance} \
        --output {params.out_folder} \
        --crisproff \
        --dict_genome_info {input.dic} \
        &>{log}
    """
