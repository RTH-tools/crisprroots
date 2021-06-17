rule RIsearch2_search:
    """
    RIsearch2 search for CRISPRoff
    """
    input:
        indexed_genome="%s/6_GATK_variants/variated_genome.suf" % config["results_folder"]
    output:
        directory("%s/6-5_CRISPRoff/RIsearch/" % config["results_folder"])
    log:
        "%s/logs/6-5-1_RIsearch2_search.log" % config["results_folder"]
    params:
        gRNA=config["Endonuclease"]["gRNA_with_PAM_fasta"],
        scripts_folder=config["path_to_snakemake"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -q : query sequence
        # -i : indexed target
        # -s : seed region
        # -m : mismatches in seed
        # -e : energy threshold
        # -l : length of sequences flanking both ends of the seed to be considered for extension
        # --noGUseed : disable wobble base pairs in seed
        # -p3 : output format
        # -t : num threads
    
        RIsearch2 -q {params.gRNA} -i {input.indexed_genome} -s 1:20 -m 6:0 -e 10000 -l 0 --noGUseed -p3 -t 16 &>{log}
        mv risearch_*.out.gz {output}
    """
