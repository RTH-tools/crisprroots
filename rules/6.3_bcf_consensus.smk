rule BCFTOOLS_consensus:
    """
    Intersect germline variant calls to get variant consensus
    """
    input:
        variants = "%s/6_GATK_variants/intersection_original/0000_filtered.vcf.gz" % config["results_folder"]
    output:
        genome = "%s/6_GATK_variants/variated_genome.fa" % config["results_folder"],
        chain_liftover = "%s/6_GATK_variants/chain_liftover_variated_genome.chain" % config["results_folder"]
    log:
        "%s/logs/6-3_bcftools_consensus.log" % config["results_folder"]
    params:
        reference=config["picard_reference"],
        heterozygous_allele_keep=config["BCF_consensus"]["heterozygous_keep"]
    conda:
        "../envs/py3.yaml"
    shell:"""
    
    #******PARAMETERS*****
    # -f : Reference fasta file
    # -o : Output file
    # -c : Write a chain liftover file
    # -H : choose which allele to keep in the consensus: 
    #       R: the REF allele (in heterozygous genotypes)
    #       A: the ALT allele (in heterozygous genotypes)

        bcftools consensus -f {params.reference} {input.variants} -o {output.genome} -c {output.chain_liftover} -H {params.heterozygous_allele_keep} &>{log}
    """
