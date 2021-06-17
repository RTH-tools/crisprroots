rule PICARD_sortaligned:
    """
    Use picard to sort the aligned reads by queryname
    """
    input:
        bam_aligned="%s/2_sortaligned/{sample}/Aligned.Sorted.bam" % config["results_folder"],
    output:
        bam_sorted=temp("%s/3_picard_sortaligned/{sample}/{sample}.Aligned.Sorted.bam" % config["results_folder"]),
    log:
        "%s/logs/3_picard_sort_aligned_{sample}.log" % config["results_folder"],
    conda:
        "../envs/gatk-picard.yaml"
    shell: """
        # sort aligned reads
        picard SortSam \
        INPUT={input.bam_aligned} \
        OUTPUT={output.bam_sorted} \
        SORT_ORDER=queryname \
        &>{log}
    """
