rule BEDTOOLS_intersect_degenes_offtargets:
    input:
        off_cutpos="%s/12-0_OffTargetCutPos/offtarget_cut_positions.bed" % config["results_folder"],
        genes_coords="%s/12-2_DeGeneCoords/genes_coordinates.bed" % config["results_folder"]
    output:
        go_intersect=temp("%s/12-3_intersect_degenes_offtargets/genes_offtargets_intersect.bed" % config["results_folder"]),
        g_sorted=temp("%s/12-3_intersect_degenes_offtargets/genes_coordinates_sorted.bed" % config["results_folder"]),
        o_sorted=temp("%s/12-3_intersect_degenes_offtargets/offtarget_cut_positions_sorted.bed" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
        out_folder=directory("%s/12-3_intersect_degenes_offtargets" % config["results_folder"])
    singularity: config["Singularity"]
    shell:
        """
        
        #******PARAMETERS*****
        # -a, -b : files to compare. Each feature in A is comapred to B. 
        # -loj : left outer join
        # -wa : write original A entries in the output
        # -sorted : files are sorted, reduces memory usage
        
        bedtools sort -i {input.off_cutpos} >{params.out_folder}/offtarget_cut_positions_sorted.bed
        bedtools sort -i {input.genes_coords} >{params.out_folder}/genes_coordinates_sorted.bed
        bedtools intersect \
        -a {params.out_folder}/offtarget_cut_positions_sorted.bed \
        -b {params.out_folder}/genes_coordinates_sorted.bed -loj -wa -sorted \
        >{params.out_folder}/genes_offtargets_intersect.bed
    """
