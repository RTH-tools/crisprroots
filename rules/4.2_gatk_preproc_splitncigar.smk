rule GATK_splitncigar:
    """
    Use GATK to split N cigar reads
    """
    input:
        dedup="%s/4_GATK_dedupSplit/{sample}/{sample}.Dedup.sorted.bam" % config["results_folder"]
    output:
        split=temp("%s/4_GATK_dedupSplit/{sample}/{sample}.Split.bam" % config["results_folder"]),
        dedup_idx=temp("%s/4_GATK_dedupSplit/{sample}/{sample}.Split.bai" % config["results_folder"])
    log:
        "%s/logs/4-2_GATK_SplitNCigar_{sample}.log" % config["results_folder"],
    params:
        reference="%s" % config["picard_reference"]
    singularity: config["Singularity"]
    shell: """
    
        #******PARAMETERS*****
        # OBI = create output bam index
        # I = input bam
        # R = reference file
        # O = output split bam 
        
        gatk SplitNCigarReads \
        -R {params.reference} \
        -I {input.dedup} \
        -O {output.split} \
        -OBI \
        &>{log}
    """
