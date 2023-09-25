rule ROOTS_collapse_coordinates:
    """
    Off-targets that overlap more than one gene/features are collapsed in one line
    """
    input:
        genes_off_intersect="%s/12-3_intersect_degenes_offtargets/genes_offtargets_intersect.bed" % config[
            "results_folder"]
    output:
        out_collapsed=temp("%s/12-4_CollapseCoordinatesGenesOff/collapsed_genes_offtargets_coordinates.bed" % config[
            "results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
    log:
        "%s/logs/12-4_CollapseCoordinatesGenesOff.log" % config["results_folder"]
    singularity: config["Singularity"]
    shell:
        """
        python3 {params.scripts_folder}/scripts/12.4_collapse_expression_offtargets.py \
        --bed_intersect {input.genes_off_intersect} \
        --out_file {output.out_collapsed} &>{log}

    """
