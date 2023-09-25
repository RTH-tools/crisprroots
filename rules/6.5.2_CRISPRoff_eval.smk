rule CRISPRoff:
    """
    Use CRISPRoff to find off-targets
    """
    input:
        "%s/6-5_CRISPRoff/RIsearch2/risearch_%s.out.gz" % (config["results_folder"],gRNA_ID)
    output:
        CRISPRspec="%s/6-5_CRISPRoff/CRISPRspec.tsv" % config["results_folder"],
        CRISPRparams="%s/6-5_CRISPRoff/CRISPRparams.tsv" % config["results_folder"],
        CRISPRoff="%s/6-5_CRISPRoff/%s.CRISPRoff.tsv" % (config["results_folder"], gRNA_with_PAM),
    log:
        "%s/logs/6-5-2_CRISPRoff_eval.log" % config["results_folder"]
    params:
        gRNA=config["gRNA_with_PAM_fasta"],
        folder="%s/6-5_CRISPRoff/" % config["results_folder"],
        scripts_folder=config["CRISPRroots"],
        gRNA_plus_PAM="%s.CRISPRoff.tsv" % gRNA_with_PAM
    singularity: config["Singularity"]
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
    
        risearch_dir="$(dirname "{input}")"
        echo $risearch_dir
        python3 \
        {params.scripts_folder}/scripts/crisproff-1.1.2/CRISPRspec_CRISPRoff_pipeline.py \
        --guides {params.gRNA} \
        --risearch_results_folder $risearch_dir \
        --no_azimuth \
        --duplex_energy_params \
        {params.scripts_folder}/scripts/crisproff-1.1.2/energy_dics.pkl \
        --specificity_report {output.CRISPRspec} \
        --guide_params_out {output.CRISPRparams} \
        --evaluate_all \
        --CRISPRoff_scores_folder {params.folder} &>{log}
        
        rm temp_grna_id_ss.ps
    """
