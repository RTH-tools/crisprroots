rule MULTIQC_compose_input:
    """
    Collect FASTQC output paths in one text file
    """
    input:
        expand("%s/preproc/0_fastqc/{sample}" % config["results_folder"],sample=lst_samples)
    output:
        txt="%s/preproc/0_fastqc/input.lst" % config["results_folder"]
    run:
        with open(output.txt,"w") as out:
            print(*input,sep="\n",file=out)
