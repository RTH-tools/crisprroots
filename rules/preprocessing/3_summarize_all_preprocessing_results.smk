rule ROOTS_summarize_stats_multiqc:
    """
    Summarize all multiqc stats
    """
    input:
        mqc_1="%s/preproc/1-2_multiqc_after_cutadapt_cleaning/multiqc_data/multiqc_general_stats.txt" % config[
            "results_folder"],
        mqc_2="%s/preproc/2-2_multiqc_after_rRNA_removal/multiqc_data/multiqc_general_stats.txt" % config[
            "results_folder"],
        mqc_0="%s/preproc/0-1_multiqc/multiqc_data/multiqc_general_stats.txt" % config["results_folder"]
    output:
        mqc_summary="%s/../report/multiqc_samples_stats.xlsx" % config["results_folder"]
    log:
        "%s/logs/preproc/3_summary_multiqc.log" % config["results_folder"]
    threads: 1
    conda:
        "../../envs/py3.yaml"
    params:
        scripts_folder=config['path_to_snakemake']
    shell: """
        printf \"Summarizing all multiqc main results\\n\"

        #******PARAMETERS*****
        # -i : input temp folder containing multiqc assessment tables
        # -o : output summary file

        output_dir=$(dirname {output.mqc_summary})
        mkdir -p $output_dir
        mkdir -p $output_dir/tmp
        cp {input.mqc_0} $output_dir/tmp/reads_initial_assessment.tsv
        cp {input.mqc_1} $output_dir/tmp/reads_retained_adapter_removed_quality_and_minLen_filtered.tsv
        cp {input.mqc_2} $output_dir/tmp/reads_retained_rna_depleted.tsv
        python3 {params.scripts_folder}/scripts/p3_gather_preproc_stats.py -i $output_dir/tmp -o {output.mqc_summary} &> {log}
        rm -r $output_dir/tmp
    """
