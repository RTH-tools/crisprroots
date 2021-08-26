rule MULTIQC_compose_input_after_cleaning:
    """
    Collect FASTQC output paths in one text file
    """
    input:
        expand("%s/preproc/1-1_fastqc_after_cutadapt_cleaning/{sample}" % config["results_folder"],sample=lst_samples)
    output:
        txt="%s/preproc/1-1_fastqc_after_cutadapt_cleaning/input.lst" % config["results_folder"]
    run:
        with open(output.txt,"w") as out:
            print(*input,sep="\n",file=out)
