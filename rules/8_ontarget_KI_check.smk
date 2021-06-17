rule ROOTS_ontarget_check_knockin:
    """
    Check for the presence of on-target mutations
    """
    input:
        bams=expand("%s/2_sortaligned/{sample}/Aligned.Sorted.bam" % config["results_folder"],sample=lst_samples),
        indexes=expand("%s/2_sortaligned/{sample}/Aligned.Sorted.bam.bai" % config["results_folder"],sample=lst_samples)
    output:
        report="%s/../report/on_target_knockin.xlsx" % config["results_folder"]
    log:
        "%s/logs/8_ROOTS_ontarget_check.log" % config["results_folder"]
    params:
        scripts_folder=config["path_to_snakemake"],
        positions=config["Edits"]["position"],
        mutants=config["Edits"]["mutant"],
        ref_fasta=config["picard_reference"],
        splice_acceptor=config["Edits"]["splice_acceptor"],
        splice_donor=config["Edits"]["splice_donor"],
        intron=config["Edits"]["intron"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -p : Genomic position(s) to analyze
        # -m : Mutation to search for, given in the same order as the corresponding positions specified in parameter -p
        # -r : Fasta file of the reference genome
        # -b : Alignment file(s) in bam format. The correspornding index should also be present
        # -o : Output file
        # -sa : The mutation affects a spice acceptor. For each mutation, use value 1 as "yes", 0 as "no" (eg. 1 0 0)
        # -sd : The mutation affects a spice donor. For each mutation, use value 1 as "yes", 0 as "no" (eg. 0 0 1)
        # -i : The mutation affects an intron. For each mutation, use value 1 as "yes", 0 as "no" (eg. 0 1 0 )

        python3 {params.scripts_folder}/scripts/8_check_KI.py \
                -p {params.positions} \
                -m {params.mutants} \
                -r {params.ref_fasta} \
                -b {input.bams} \
                -o {output.report} \
                -sa {params.splice_acceptor} \
                -sd {params.splice_donor} \
                -i {params.intron} \
                &>{log}
    """
