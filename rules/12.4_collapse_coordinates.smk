rule ROOTS_collapse_coordinates:
    """
    Off-targets that overlap more than one gene/features are collapsed in one line
    """
    input:
        genes_off_intersect="%s/12-3_intersect_degenes_offtargets/genes_offtargets_intersect.bed" % config[
            "results_folder"]
    output:
        out_collapsed="%s/12-4_CollapseCoordinatesGenesOff/collapsed_genes_offtargets_coordinates.bed" % config[
            "results_folder"]
    params:
        scripts_folder=config["path_to_snakemake"],
    log:
        "%s/logs/12-4_CollapseCoordinatesGenesOff.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell:
        """
        python3 {params.scripts_folder}/scripts/12.4_collapse_expression_offtargets.py \
        --bed_intersect {input.genes_off_intersect} \
        --out_file {output.out_collapsed} &>{log}

    """
