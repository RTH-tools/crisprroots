rule ROOTS_dfgenecoordinates:
    """
    Get genes coordinates and expand to promoters
    """
    input:
        dic="%s/0_utils/dict_chroms_lengths.pkl" % config["results_folder"],
        lifted_annotations="%s/0_utils/lifted_annotations.gtf" % config["results_folder"],
        diffexp_table="%s/12-1_DESEq2DE/diffexpr-results.tsv" % config["results_folder"]
    output:
        "%s/12-2_DeGeneCoords/genes_coordinates.bed" % config["results_folder"]
    params:
        dir="%s/12-1_DESEq2DE" % config["results_folder"],
        scripts_folder=config["path_to_snakemake"],
        len_promoter=config["ExpressionBasedScreening"]["len_promoter"]
    log:
        "%s/logs/12-2_degenescoordinates.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell:
        """
        #******PARAMETERS*****
        # -annotations : Annotations file in gtf format
        # -diffexp_genes : Table containing strongly differentially expressed genes to be analyzed
        # -o : Output file name
        # --length : Length of the promoter region
        # -d : Path to dictionary of chromosome length
        
        arr=($(ls {params.dir}/genes_*.tsv))

        python3 {params.scripts_folder}/scripts/12.2_get_gene_coords.py \
        --annotations {input.lifted_annotations} \
        --diffexp_genes "${{arr[@]}}" \
        -o {output} \
        --length {params.len_promoter} \
        -d {input.dic} \
        &>{log}

    """
