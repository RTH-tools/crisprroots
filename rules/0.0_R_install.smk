rule R_install:
    output:
        Rlog=temp("%s/0_utils/0-0_R_install_check.tmp" % config["results_folder"])
    log:
        R_install="%s/logs/0-0_R_install.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    params:
        scripts_folder=config["path_to_snakemake"],
        threads=config["R_install_pkgs"]["threads"]
    shell: """

        #******PARAMETERS*****
        # (1) : num threads
        # (2) : Rlog check output file

        Rscript {params.scripts_folder}/scripts/0.0_install_R_pkgs.R {params.threads} {output.Rlog} &>{log.R_install}
"""
