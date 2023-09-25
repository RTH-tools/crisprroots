rule ROOTS_getOfftargetCutPositions:
    """
    Get the coordinates of the possible cut positions from the CRISPRoff output
    OR from a bed file with predicted off-targets. CRISPRoff table version.
    """
    input:
        offtargets="%s/6-5_CRISPRoff/%s.CRISPRoff.tsv" % (config["results_folder"], gRNA_with_PAM)
    output:
        off_cutpos=temp("%s/12-0_OffTargetCutPos/offtarget_cut_positions.bed" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
        cp=config["Endonuclease"]["cut_position"],
        gRNA=config["Endonuclease"]["gRNA_sequence"],
        report=config["report_folder"]
    log:
        "%s/logs/12-0_Off_target_cut_pos.log" % config["results_folder"]
    singularity: config["Singularity"]
    shell: """
    
        #******PARAMETERS*****
        # -i : Path to input file 
        # -o : Path to output file
        # -cp : num. nts from the 3' of the predicted off-target coordinates to the cleavage position
        # --crisproff : Flag. Use if the input is the crisproff tsv output table
        # --webserver : Flag. Use if the input is from the CRISPRoff webserver
        # -g : gRNA sequence
        # -g : Path to report folder, used for alerts (eg. two on-target)
    
        python3 {params.scripts_folder}/scripts/12.0_get_offtarget_cutpos.py \
        --input {input.offtargets} \
        --cut_position {params.cp} \
        --crisproff \
        --output {output.off_cutpos} \
        -g {params.gRNA} \
        -f {params.report} \
        &>{log}
    """
