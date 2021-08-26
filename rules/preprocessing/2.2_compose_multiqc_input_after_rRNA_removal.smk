rule MULTIQC_compose_input_after_rRNA_removal:
    """
    Collect FASTQC output paths in one text file
    """
    input:
        expand("%s/preproc/2-1_fastqc_after_rRNA_removal/{sample}" % config["results_folder"],sample=lst_samples)
    output:
        txt="%s/preproc/2-1_fastqc_after_rRNA_removal/input.lst" % config["results_folder"]
    run:
        with open(output.txt,"w") as out:
            print(*input,sep="\n",file=out)
