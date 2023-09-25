rule GATK_mutect2_chromosome_wise:
    """
    Run MuTect2 chromosome-wise parallelly (using GNU parallel)
    MuTect2 calls somatic variants, allowing for multiple different allelic fractions for each variant
    """
    input:
        edited=expand("%s/5_chromosome_wise_splitted_bam/{edited}/" % config["results_folder"],edited=lst_edited),
        original=expand("%s/5_chromosome_wise_splitted_bam/{original}/" % config["results_folder"],original=lst_original),
        chrsplit_checkpoint1=expand("%s/5_chromosome_wise_splitted_bam/{edited}/split.done" % config["results_folder"],edited=lst_edited),
        chrsplit_checkpoint2=expand("%s/5_chromosome_wise_splitted_bam/{original}/split.done" % config["results_folder"],original=lst_original)
    output:
        checkpoint="%s/7_GATK_mutect2_chromosome_wise/mutect2.done" % config["results_folder"],
    log:
        "%s/logs/7_GATK_mutect2_chromosome_wise.log" % config["results_folder"],
    params:
        outdir="%s/7_GATK_mutect2_chromosome_wise/" % config["results_folder"],
        indir="%s/5_chromosome_wise_splitted_bam/" % config["results_folder"],
        BQST=config["Mutect2"]["base_quality_score_threshold"],
        MQB=config["Mutect2"]["min_base_quality_score"],
        CD=config["Mutect2"]["callable_depth"],
        scripts_folder=config["CRISPRroots"],
        reference=config["picard_reference"],
        parallel_jobs=config["Mutect2"]["num_threads"]

    singularity: config["Singularity"]
    shell: """
    
        #******PARAMETERS*****
        # Mutect2 params in script 7_run_gatk.py
        # -bqst : --base qualiy score threshold. Base qualities below this threshold will be reduced to the minimum (6)
        # -cd : --callable depth. Minimum depth to be considered callable for Mutect stats. Does not affect genotyping.
        # -mqb : --min-base-quality-score: Minimum base quality required to consider a base for calling,
        # -r : path to reference genome
        # -ed : path to list of edited samples
        # -or : path to list of original samples 
        # -O : Path to filtered vcf output
        # -o : Path to vcf output

        infile_edited=$(for i in {input.edited}; do echo $(find $i -type f -name '*.bam'); done)
        infile_original=$(for i in {input.original}; do echo $(find $i -type f -name '*.bam'); done)
        tot_samples=$(find {params.indir}/* -type d | wc -l)
        echo ""
        echo "Total number of samples (edited/original): $tot_samples"
        infile_edited_chr=$({params.scripts_folder}/scripts/7_find_common_chr_btw_samples.sh {params.indir})
        
        echo ""      
        echo MuTect2 will run for the following : $infile_edited_chr 

        # function to invoke script calling gatk mutect2 parallelly 
        function invoke_mutect2_script {{

           x=$1

           infile_edited=$(for i in {input.edited}; do echo $(find $i -type f -name '*.bam' | grep $x); done);
           infile_original=$(for i in {input.original}; do echo $(find $i -type f -name '*.bam' | grep $x); done);

           out=$(echo {params.outdir}$x".vcf" | sed 's/\.bam//g;')
           out2=$(echo {params.outdir}$x"_filtered.vcf" | sed 's/\.bam//g;')
           logout=$(echo {log} | sed "s/\.log/_$x/g; s/.bam/.log/g")

           echo \ "output file name for $x:" $out 
           echo \ input edited: $infile_edited 
           echo \ input original: $infile_original 
           echo "Running GATK_mutect2..."
           echo ""                               
            
           python3 {params.scripts_folder}/scripts/7_run_gatk.py \
           -r {params.reference} -bqst {params.BQST} \
           -cd {params.CD} -mqb {params.MQB} -o $out \
           -ed $infile_edited -or $infile_original -O $out2 2>$logout

         }}

        # calling GNU parallel
        export -f invoke_mutect2_script
        parallel --verbose -j {params.parallel_jobs} invoke_mutect2_script ::: $infile_edited_chr

        touch "{output.checkpoint}"

   """
