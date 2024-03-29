rule BEDTOOLS_intersect_variants_genes:
    """
    Use bedtools to intesect the list of variants called with genes annotations (lifted)
    """
    input:
        variants="%s/11_VariantBasedScreening/EvaluatedVariantsOffTargets.bed" % config["results_folder"],
        lifted_annotations="%s/0_utils/lifted_annotations.gtf" % config["results_folder"]
    output:
        v_g_intersect=temp("%s/11_VariantBasedScreening/variants_genes_intersect.bed" % config["results_folder"]),
        g_sort=temp("%s/11_VariantBasedScreening/lifted_genes_coordinates_sorted.bed" % config["results_folder"]),
        v_sort=temp("%s/11_VariantBasedScreening/variants_cut_positions_sorted.bed" % config["results_folder"])
    params:
        scripts_folder=config["CRISPRroots"],
        out_folder=directory("%s/11_VariantBasedScreening" % config["results_folder"])
    singularity: config["Singularity"]
    shell:
        """
        
        #******PARAMETERS*****
        # -a, -b : files to compare. Each feature in A is comapred to B. 
        # -loj : left outer join
        # -wa : write original A entries in the output
        # -sorted : files are sorted, reduces memory usage
        
        bedtools sort -i {input.variants} >{params.out_folder}/variants_cut_positions_sorted.bed
        bedtools sort -i {input.lifted_annotations} >{params.out_folder}/lifted_genes_coordinates_sorted.bed
        bedtools intersect \
        -a {params.out_folder}/variants_cut_positions_sorted.bed \
        -b {params.out_folder}/lifted_genes_coordinates_sorted.bed -loj -wa -sorted \
        >{output.v_g_intersect}
    """
