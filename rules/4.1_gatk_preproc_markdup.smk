rule GATK_markdup:
    """
    Use GATK to mark duplicates (and sort)
    """
    input:
        bam_merged="%s/3_picard_sortaligned/{sample}/{sample}.Aligned.Sorted.bam" % config["results_folder"],
    output:
        dedup=temp("%s/4_GATK_dedupSplit/{sample}/{sample}.Dedup.bam" % config["results_folder"]),
        dedup_sorted=temp("%s/4_GATK_dedupSplit/{sample}/{sample}.Dedup.sorted.bam" % config["results_folder"]),
        metrics=temp("%s/4_GATK_dedupSplit/{sample}/{sample}.Dedup.Metrics" % config["results_folder"]),
    log:
        "%s/logs/4-1_GATK_MarkDup_{sample}.log" % config["results_folder"],
    singularity: config["Singularity"]
    shell: """
        #marking duplicate reads

        gatk MarkDuplicates \
        --INPUT {input.bam_merged} \
        --OUTPUT {output.dedup} \
        --METRICS_FILE {output.metrics} \
        &>{log}

        picard SortSam \
        INPUT={output.dedup} \
        OUTPUT={output.dedup_sorted} \
        SORT_ORDER=coordinate
    """
