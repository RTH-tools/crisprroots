rule SAMTOOLS_sort_index_bam:
    input:
        bam_aligned="%s/1_star_align2pass/{sample}/{sample}.Aligned.out.bam" %config["results_folder"]
    output:
        sorted_bai="%s/2_sortaligned/{sample}/Aligned.Sorted.bam.bai" %config["results_folder"],
        sorted_bam="%s/2_sortaligned/{sample}/Aligned.Sorted.bam" %config["results_folder"],
    log:
        sort="%s/logs/2_sort_bam_{sample}.log" % config["results_folder"],
        index="%s/logs/2_index_bam_{sample}.log" % config["results_folder"],
    conda:
         "../envs/py3.yaml"
    shell:"""
        samtools sort {input.bam_aligned} -o {output.sorted_bam} &>{log.sort}
        samtools index {output.sorted_bam} &>{log.index}

    """