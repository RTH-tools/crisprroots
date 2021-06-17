rule RIsearch2_indexing:
    """
    RIsearch2 index the variated genome
    """
    input:
        genome="%s/6_GATK_variants/variated_genome.fa" % config["results_folder"]
    output:
        indexed_genome=temp("%s/6_GATK_variants/variated_genome.suf" % config["results_folder"])
    log:
        "%s/logs/6-4_RIsearch2_index.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
        
        #******PARAMETERS*****
        # -c : path to target fasta sequence to be indexed
        # -o : path to output file
        
        RIsearch2 -c {input.genome} -o {output.indexed_genome} &>{log}
    """
