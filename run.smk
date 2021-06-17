import pandas as pd
from snakemake.utils import validate
import os

#***** load config and sample sheets *****#
configfile: "./config.yaml"
validate(config,schema="schemas/config.schema.yaml")

df_samples = pd.read_csv(config["samples_table"],sep="\t",index_col="Sample_ID")
invalid_conditions = [x for x in df_samples["Condition"].unique() if x not in ["Original", "Edited"]]
assert len(invalid_conditions) == 0, "The samples table contains invalid conditions. Only conditions \"Edited\" and \"Original\" are allowed."

lst_edited = list(df_samples[df_samples['Condition'] == 'Edited'].index)
lst_original = list(df_samples[df_samples['Condition'] == 'Original'].index)
lst_samples = list(df_samples.index)
with open(config["Endonuclease"]["gRNA_with_PAM_fasta"],'r') as gRNA_file:
    gRNA_with_PAM = gRNA_file.readlines()[1].rstrip('\n')

if not os.path.exists(config["results_folder"]):
    os.mkdir(config["results_folder"])
if not os.path.exists("%s/0_utils" % config["results_folder"]):
    os.mkdir("%s/0_utils" % config["results_folder"])
if not os.path.exists("%s/0_utils/edits.bed" % config["results_folder"]):
    with open("%s/0_utils/edits.bed" % config["results_folder"],"w") as out:
        for edit in config["Edits"]["position"]:
            chr, pos = edit.split(':')
            out.write("%s\t%i\t%s\n" % (chr, int(pos) - 1, pos))

#***** run pipeline *****#
rule all:
    input:
        "%s/../report/candidate_off_targets.xlsx" % config["results_folder"],
        "%s/../report/multiqc_samples_stats.xlsx" % config["results_folder"],
        "%s/../report/mapping_stats.xlsx" % config["results_folder"],
        "%s/../report/on_target_knockin.xlsx" % config["results_folder"],
        "%s/../report/on_target_knockout.xlsx" % config["results_folder"],
        expand("%s/2-2_RSeQC_libtype/{sample}_libtype.txt" % config["results_folder"],sample=lst_samples),
    #expand("%s/eSNPKaryotyping/{sample}_Zygosity_Blocks.pdf" % config["report_folder"], sample=lst_samples)

rule off_targets:
    input:
        txt="%s/../report/candidate_off_targets.xlsx" % config["results_folder"]

rule preproc_and_map:
    input:
        "%s/../report/multiqc_samples_stats.xlsx" % config["results_folder"],
        "%s/../report/mapping_stats.xlsx" % config["results_folder"],
        expand("%s/2_sortaligned/{sample}/Aligned.Sorted.bam" % config["results_folder"],sample=lst_samples)

rule variants_to_genome:
    input:
        expand("%s/6_GATK_variants/{sample}/variants_filtered.vcf" % config["results_folder"],sample=lst_samples)

rule eSNPKaryotyping:
    input:
        expand("%s/eSNPKaryotyping/{sample}_Zygosity_Blocks.pdf" % config["report_folder"],sample=lst_samples)

rule on_target_check:
    input:
        "%s/../report/on_target_knockin.xlsx" % config["results_folder"],
        "%s/../report/on_target_knockout.xlsx" % config["results_folder"]

rule get_lib_type:
    input:
        expand("%s/2-2_RSeQC_libtype/{sample}_libtype.txt" % config["results_folder"],sample=lst_samples)

rule get_variated_genome:
    input:
        "%s/6_GATK_variants/variated_genome.fa" % config["results_folder"]


#***** load rules *****#
if config["sequencing"] == "single":
    include: "rules/1_star_align2pass_single.smk"
else:
    include: "rules/1_star_align2pass_paired.smk"
include: "rules/0.0_R_install.smk"
include: "rules/0.1_make_utils.smk"
include: "rules/0.2_RNAfold.smk"
include: "rules/1.1_summarize_mapping_stats.smk"
include: "rules/2_sort_index_bam.smk"
include: "rules/2.2_RSeQC_libtype.smk"
include: "rules/3_picard_sortaligned.smk"
include: "rules/4.1_gatk_preproc_markdup.smk"
include: "rules/4.2_gatk_preproc_splitncigar.smk"
include: "rules/5_split_bam_chromosome_wise.smk"
include: "rules/6_gatk_haplotypecaller.smk"
include: "rules/6.1_gatk_filter.smk"
if config["variated_genome"] == 'no':
    include: "rules/6.2_bcf_novariants.smk"
    include: "rules/12.0_offtargets_cutPos_novariants.smk"
else:
    include: "rules/6.2_bcf_intersect_original.smk"
    include: "rules/12.0_offtargets_cutPos.smk"
include: "rules/6.3_bcf_consensus.smk"
include: "rules/6.4_RIsearch2_indexing.smk"
include: "rules/6.5.1_RIsearch2_search.smk"
include: "rules/6.5.2_CRISPRoff_eval.smk"
include: "rules/7_gatk_mutect2_chromosome_wise.smk"
include: "rules/7.1_gatk_mutect2_chromosome_wise_vcf_concat.smk"
include: "rules/8_ontarget_KI_check.smk"
include: "rules/8_ontarget_KO_check.smk"
include: "rules/9_featurecounts_quantification.smk"
include: "rules/10_bedops_vcftobed.smk"
include: "rules/10.1_bedops_liftover.smk"
include: "rules/11.0_variant_screening.smk"
include: "rules/11.1_intersect_variants_genes.smk"
include: "rules/12.1_DESeq2_diffexp.smk"
include: "rules/12.2_genes_coordinates.smk"
include: "rules/12.3_genes_offTargets_intersection.smk"
include: "rules/12.4_collapse_coordinates.smk"
include: "rules/12.5_singlechr_filter.smk"
include: "rules/12.6_expression_screening.smk"
include: "rules/13_flags.smk"
include: "rules/14_eSNPKAryotyping.smk"
include: "rules/15_offTargets_report.smk"
#***** load preprocessing *****#
include: "rules/preprocessing/%s/1_cleaning_cutadapt.smk" % config["sequencing"]
include: "rules/preprocessing/%s/2_filter_rrna.smk" % config["sequencing"]
include: "rules/preprocessing/0.1_compose_multiqc_input.smk"
include: "rules/preprocessing/0.1_multiqc.smk"
include: "rules/preprocessing/%s/0_qualitycheck_fastqc.smk" % config["sequencing"]
include: "rules/preprocessing/%s/1.1_qualitycheck_after_cleaning.smk" % config["sequencing"]
include: "rules/preprocessing/1.2_compose_multiqc_input_after_cleaning.smk"
include: "rules/preprocessing/1.2_multiqc_after_cleaning.smk"
include: "rules/preprocessing/%s/2.1_qualitycheck_after_rRNA_removal.smk" % config["sequencing"]
include: "rules/preprocessing/2.2_compose_multiqc_input_after_rRNA_removal.smk"
include: "rules/preprocessing/2.2_multiqc_after_rRNA_removal.smk"
include: "rules/preprocessing/3_summarize_all_preprocessing_results.smk"
