rule CRISPRoff:
    """
    Use CRISPRoff to find off-targets
    """
    input:
        "%s/6-5_CRISPRoff/RIsearch/" % config["results_folder"]
    output:
        CRISPRspec="%s/6-5_CRISPRoff/CRISPRspec.tsv" % config["results_folder"],
        CRISPRparams="%s/6-5_CRISPRoff/CRISPRparams.tsv" % config["results_folder"],
        CRISPRoff="%s/6-5_CRISPRoff/%s.CRISPRoff.tsv" % (config["results_folder"], gRNA_with_PAM),
    log:
        "%s/logs/6-5-2_CRISPRoff_eval.log" % config["results_folder"]
    params:
        gRNA=config["Endonuclease"]["gRNA_with_PAM_fasta"],
        folder="%s/6-5_CRISPRoff/" % config["results_folder"],
        scripts_folder=config["path_to_snakemake"],
        gRNA_plus_PAM="%s.CRISPRoff.tsv" % gRNA_with_PAM
    conda:
        "../envs/py2.yaml"
    shell: """
    
        #******PARAMETERS*****
        # --guides : path to gRNA in fasta format
        # --risearch_results_folder : path to folder where to look for the research results
        # --no_azimuth : do not produce azimuth scores
        # --duplex_energy_parameters : path to the file with the duplex energy parameters
        # --specificity_report : output path for the gRNA specificity report
        # --guide_params_out : output path for the gRNA energy parameters
        # --evaluate_all : evaluate the gRNA even if its on-target is not in the genome
        # --CRISPRoff_scores_folder : path where the CRISPRoff results will be stored
    
        python2 \
        {params.scripts_folder}/scripts/crisproff-1.1.1/CRISPRspec_CRISPRoff_pipeline.py \
        --guides {params.gRNA} \
        --risearch_results_folder {input} \
        --no_azimuth \
        --duplex_energy_params \
        {params.scripts_folder}/scripts/crisproff-1.1.1/energy_dics.pkl \
        --specificity_report {output.CRISPRspec} \
        --guide_params_out {output.CRISPRparams} \
        --evaluate_all \
        --CRISPRoff_scores_folder {params.folder} &>{log}
        
        rm temp_grna_id_ss.ps
    """
