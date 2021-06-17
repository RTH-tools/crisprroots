rule ROOTS_singlechr_filter:
    """
    Filters off-targets in single-allele genes, after checking that no called variant overlaps the region
    """
    input:
        collapsed_coords="%s/12-4_CollapseCoordinatesGenesOff/collapsed_genes_offtargets_coordinates.bed" % config[
            "results_folder"],
    output:
        collapsed_coords="%s/12-4_CollapseCoordinatesGenesOff/collapsed_genes_offtargets_coordinates_singlechr_filtered.bed" %
                         config["results_folder"]
    params:
        scripts_folder=config["path_to_snakemake"],
        single_chr=config["single_chr"]
    log:
        "%s/logs/12-4_SingleChr_filter.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -c : Path to collapsed off-target coordinates 
        # -o : Path to output file
        # -s : List of single chromosomes (eg. chrY)
        
        python3 {params.scripts_folder}/scripts/12.5_filter_single_chr.py \
        -c {input.collapsed_coords} \
        -o {output.collapsed_coords} \
        -s {params.single_chr} \
        &>{log}
    """
