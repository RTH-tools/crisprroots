rule ROOTS_offTargetReport:
    """
    Merge expression-based and variant-based off-target reports into one report, apply filters, categorization, sorting
    """
    input:
        ofv="%s/11_VariantBasedScreening/EvaluatedVariantsOffTargets.tsv" % config["results_folder"],
        vars_genes="%s/11_VariantBasedScreening/variants_genes_intersect.bed" % config["results_folder"],
        oft="%s/12-5_ExpressionBasedScreening/EvaluatedExpressionOffTargets.tsv" % config["results_folder"],
        frt="%s/13_flags/frt.bed" % config["results_folder"],
        frv="%s/13_flags/frv.bed" % config["results_folder"],
        fsv="%s/13_flags/fsv.bed" % config["results_folder"]
    output:
        txt="%s/../report/candidate_off_targets.xlsx" % config["results_folder"]
    log:
        "%s/logs/15_Offtarget_report.log" % config["results_folder"]
    params:
        scripts_folder=config["CRISPRroots"],
        seed=config["Endonuclease"]["seed_region"],
        g=config["Endonuclease"]["gRNA_sequence"],
        cp=config["Endonuclease"]["cut_position"],
        extend_binding=config["Endonuclease"]["extend_binding"],
        gRNA=config["Endonuclease"]["gRNA_sequence"],
        max_seed_mismatch=config["Endonuclease"]["max_mm_seed"],
        deltagb=config["Endonuclease"]["eng_threshold"],
        targets=config["Edits"]["KO"],
        type=config["Edits"]["type"],
        resource=config["CRISPRroots"]+'/resources/candidate_off_targets_template.xlsx'
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -oft : table with potential off-targets from expression-based off-target screening
        # -ofv : table with potential off-targets from variant-based off-target screening
        # -seed : length of the seed region
        # -gRNA_sequence : gRNA sequence
        # --extend_binding : max size allowed for extend_binding
        # --cut_position : Number of nucleotides from the 3p end of the predicted off-target coordinates to the cleavage position
        # --seed_mismatch_tolerance : Max mismatches allowed in seed
        # --float_max_deltagb : Maximum weighted binding energy allowed in gRNA-DNA binding
        # --flag_repeatmask_transcripts : Path to bed file with off-targets from expression-based analysis flagged because repeatmasked
        # --flag_repeatmask_variants : Path to bed file with off-targets from variant-based analysis flagged because repeatmasked
        # --flag_dbsnp_variants : Path to bed file with off-targets from variant-based analysis flagged because overlapping a dbSNP entry
        # --vars_genes : bed file of coordinates of genes intersecting variants
        # -o : Path to output file
        # -ko : Knockout gene IDs
        # -type : Type of edit: KI (knockin) or KO (knockout)
   
        python3 {params.scripts_folder}/scripts/15_generate_report.py \
        -oft {input.oft} \
        -ofv {input.ofv} \
        --seed {params.seed} \
        --cut_position {params.cp} \
        --gRNA_sequence {params.g} \
        --extend_binding {params.extend_binding} \
        --seed_mismatch_tolerance {params.max_seed_mismatch} \
        --float_max_deltagb {params.deltagb} \
        --flag_repeatmask_transcripts {input.frt} \
        --flag_repeatmask_variants {input.frv} \
        --flag_dbsnp_variants {input.fsv} \
        --vars_genes {input.vars_genes} \
        -type {params.type} \
        -ko {params.targets} \
        -res {params.resource} \
        -o {output} &>{log}
        
    """
