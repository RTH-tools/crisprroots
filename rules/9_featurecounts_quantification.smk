rule FEATURECOUNTS_quantification:
    """
    Quantify gene expression (FeatureCounts)
    """
    input:
        bam_aligned=expand("%s/1_star_align2pass/{sample}/{sample}.Aligned.out.bam" % config[
            "results_folder"],sample=lst_samples)
    output:
        featurecounte_out="%s/9_featurecounts_quantification/feature_count.tsv" % config["results_folder"],
        featurecounts_stats="%s/9_featurecounts_quantification/feature_count.tsv.summary" % config["results_folder"]
    params:
        libtype=config["Featurecounts"]["libtype"],
        seq_type=config["sequencing"],
        annot=config["annotations_gtf"]
    log:
        "%s/logs/9_featurecounts_quantification.log" % config["results_folder"]
    threads: 16
    conda:
        "../envs/star-featurecounts.yaml"
    shell: """
        printf \"Feature counts\\n\"

        #******PARAMETERS*****
        #-p = paired end
        #-C = do not count chimeric fragments (paired fragments having ends aligned to different chromosomes
        #-B = Only count reads with both ends mapped
        #-s = specify strand: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default 0. For paired-end reads, strand of the first read is taken as the strand of the whole fragment.

        PAIRED_OPTION=""
        if [ {params.seq_type} == "paired" ]
        then
            PAIRED_OPTION="-p"
            echo The paired option $PAIRED_OPTION was selected
        fi
        featureCounts -T {threads} -F GTF \
        -s {params.libtype} \
        $PAIRED_OPTION -C -B --extraAttributes gene_name \
        -a {params.annot} -o {output.featurecounte_out} \
        {input.bam_aligned} &> {log}
    """
