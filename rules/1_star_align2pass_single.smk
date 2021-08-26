rule STAR_align2pass:
    """
    Align reads to genome with STAR
    """
    input:
        reads=[
            "%s/preproc/2_bbduk_rrna_filter/{sample}/{sample}%s" % (config["results_folder"], config["sample_suffix"])],
    output:
        bam_aligned=temp("%s/1_star_align2pass/{sample}/{sample}.Aligned.out.bam" % config["results_folder"]),
        stats=temp("%s/1_star_align2pass/{sample}/{sample}.Log.final.out" % config["results_folder"])
    log:
        "%s/logs/1_star_align2pass_{sample}.log" % config["results_folder"]
    params:
        threads=config["STAR"]["threads"],
        dir_prefix=directory("%s/1_star_align2pass/{sample}" % config["results_folder"]),
        indexed_transcriptome=config["STAR_indexed_transcriptome"]
    conda:
        "../envs/star-featurecounts.yaml"
    shell: """
        #******PARAMETERS*****
        # --readFilesCommand : how to open the zipped format zcat
        # --genomeDir : indexed genome dir
        # --readFilesIn : input reads file
        # --outFileNamePrefix : prefix of files; includes the out folder
        # --runThreadN :  num threads to use
        # --outSAMtype : output type
        # --outSAMattrRGline : add read group info

        STAR \
        --twopassMode 'Basic' \
        --genomeDir {params.indexed_transcriptome} \
        --readFilesIn {input.reads} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.dir_prefix}/{wildcards.sample}. \
        --runThreadN {params.threads} \
        --outSAMtype BAM Unsorted \
        --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
        &>{log}
    """
