rule RSEQC_libtype:
    """
    Discovery the library type with RSEQC
    """
    input:
        bam="%s/2_sortaligned/{sample}/Aligned.Sorted.bam" % config["results_folder"]
    output:
        libtype="%s/2-2_RSeQC_libtype/{sample}_libtype.txt" % config["results_folder"]
    params:
        gene_model=config["RSeQC_gene_model"]
    conda:
        "../envs/preproc-qc.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -r : reference gene model
        # -i : bam input
    
        infer_experiment.py -r {params.gene_model} -i {input.bam} &>{output.libtype}
    """
