rule RNAfold_get_eng:
    """
    Creates the RNAfold output for the gRNA
    """
    input:
        genome="%s/6_GATK_variants/variated_genome.fa" % config["results_folder"],
    output:
        gRNA_fold="%s/0_utils/guide_RNA.fold" % config["results_folder"],
        gRNA_fa="%s/0_utils/gRNA.fa" % config["results_folder"],
    log:
        rnafold="%s/logs/0_RNAfold.log" % config["results_folder"],
    conda:
        "../envs/py2.yaml"
    params:
        gRNA=config["Endonuclease"]["gRNA_sequence"],
        scripts_folder=config["path_to_snakemake"],
    shell: """
        
        #******PARAMETERS*****
        # of : path to RNAfold output file
        # og : path to gRNA fasta file
        # g : gRNA sequence
    
        python2 {params.scripts_folder}/scripts/0.2_get_RNAfold.py \
        -og {output.gRNA_fa} \
        -of {output.gRNA_fold} \
        -g {params.gRNA} \
        &>{log.rnafold}

        rm gRNA_ss.ps
"""
