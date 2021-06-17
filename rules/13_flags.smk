rule FLAGS:
    """
    Add flags for repeatmasked and known SNPs
    """
    input:
        bed_v="%s/11_VariantBasedScreening/EvaluatedVariantsOffTargets.bed" % config["results_folder"],
        bed_t="%s/12-5_ExpressionBasedScreening/EvaluatedExpressionOffTargets.bed" % config["results_folder"],
        lifted_snpdb="%s/0_utils/lifted_SNPdb.bed" % config["results_folder"],
        lifted_rpkm="%s/0_utils/lifted_repeatmask.bed" % config["results_folder"],
    output:
        fsv="%s/13_flags/fsv.bed" % config["results_folder"],
        frv="%s/13_flags/frv.bed" % config["results_folder"],
        frt="%s/13_flags/frt.bed" % config["results_folder"]
    log:
        fsv="%s/logs/13_flags_fsv.log" % config["results_folder"],
        frt="%s/logs/13_flags_frt.log" % config["results_folder"],
        frv="%s/logs/13_flags_frv.log" % config["results_folder"]
    conda:
        "../envs/py3.yaml"
    shell: """
    
        #******PARAMETERS*****
        # -u : write A feature once if overlaps any in B
        # -a, -b : files to compare. Each feature in A is comapred to B 
        
        intersectBed \
        -u \
        -a {input.bed_v} \
        -b {input.lifted_snpdb} \
        1>{output.fsv} \
        2>{log.fsv} || touch {output.fsv}

        intersectBed \
        -u \
        -a {input.bed_v} \
        -b {input.lifted_rpkm} \
        1>{output.frv} \
        2>{log.frv} || touch {output.frv}

        intersectBed \
        -u \
        -a {input.bed_t} \
        -b {input.lifted_rpkm} \
        1>{output.frt} \
        2>{log.frt} || touch {output.frt}
    """
