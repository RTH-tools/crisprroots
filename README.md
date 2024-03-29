# CRISPRroots

CRISPRroots v.1.3: **CRISPR**–Cas9-mediated
edits with accompanying **R**NA-seq data assessed for **o**n-target and **o**ff-**t**arget **s**ites

The CRISPR/Cas9 genome editing tool can be used to study genomic variants and gene knockouts.
By combining CRISPR/Cas9 mediated editing with transcriptomic analyses it is possible to measure
the effects of genome alterations on gene expression. In such experiments it is crucial to understand
not only if the editing was successful but also if the observed differential gene expression is the result
of intended genome edits and not that of unwanted off-target effects. Potential off-targets sites that
CRISPR/Cas9 may have edited need therefore to be verified by sequencing. However, a CRISPR/Cas9
gRNA may have hundreds (or thousands) potential off-target binding sites in a genome. How to decide
which of them should be validated with higher priority? The RNA-seq data that can be sequenced
as part of a CRISPR/Cas9 experiment contains information about the sequence and expression level
of potential off-target sites/genes located in transcribed regions. Here we preset CRISPRroots, a
method that combines CRISPR/Cas9 and guide RNA binding properties, gene expression changes,
and sequence variants between edited and non-edited cells, to discover and rank potential off-targets.
The method is described in the corresponding publication (see below).

## Installation and configuration

### Prerequisites
* Snakemake ≥ 7.28.3
* Apptainer ≥ 1.1.8
* CRISPRroots docker container


The CRISPRroots pipeline can be executed via 
[Snakemake](https://snakemake.readthedocs.io/en/stable/).
All other software requirements are satisfied by the CRISPRroots container, available at [Docker Hub](https://hub.docker.com/r/gcorsi1993/crisprroots):
```shell
docker pull gcorsi1993/crisprroots
```
The container is provided to Snakemake using the option `--use-singularity`.

The pipeline was tested in a x86 64 GNU/Linux environment with Ubuntu v.18.04.1 (or newer) 
installed.

### Setup

CRISPRroots is available to download from 
[https://rth.dk/resources/crispr/](https://rth.dk/resources/crispr/). 
A test dataset is also provided on the same website. 
After downloading and un-packing the software and the test dataset, we recommend 
to explore the directory structure of the test dataset.

*  CRISPRroots_test_dataset
    *   `resources` # folder with reference files
    *   `QPRT_DEL268T_chr16_10M-40M` # folder with data examples
    *   `make_config.py` # setup script

the script `make_config.py` can be used to automatically set the paths to the data, resources, code directory and singularity in the configuration file `config.yaml` located in 
`CRISPRroots_test_dataset/QPRT_DEL268T_chr16_10M-40M`. 
To run the script, you need python3.
To setup the config file for the test dataset, run:
```shell
  $ cd CRISPRroots_test_dataset
  $ python3 make_config.py --CRISPRroots <path_to_CRISPRroots> --singularity <path_to_singularity>
```
The config file contains the parameters defined for the execution of the various steps of the pipeline. 
A copy of the config file for the CRISPRroots test dataset is provided together with 
the CRISPRroots software package, under `resources`. This file can be used as template to create the configuration file 
for your own dataset.

For a complete list of parameters, options, and usage examples for CRISPRroots please read the 
*CRISPRroots_Manual.pdf* included in the CRISPRroots software folder.

## Pipeline usage
### Basic Usage
You can run the pipeline from within the directory containing the config file (in the test dataset this is
the subfolder *QPRT_DEL268T_chr16_10M-40M*) as follows: 
```shell
$ cd QPRT_DEL268T_chr16_10M-40M # example with test dataset
$ snakemake -s <path_to_CRISPRroots>/run.smk --cores <int> --use-singularity --singularity-args
"--bind <global_path_to_QPRT_DEL268T_chr16_10M-40M>:<global_path_to_QPRT_DEL268T_chr16_10M-40M>
--bind <global_path_to_resources>:<global_path_to_resources>
--bind <global_path_to_CRISPRroots-1.3>:<global_path_to_CRISPRroots-1.3>" --dry-run
$ snakemake -s <path_to_CRISPRroots>/run.smk --cores <int> --use-singularity --singularity-args
"--bind <global_path_to_QPRT_DEL268T_chr16_10M-40M>:<global_path_to_QPRT_DEL268T_chr16_10M-40M>
--bind <global_path_to_resources>:<global_path_to_resources>
--bind <global_path_to_CRISPRroots-1.3>:<global_path_to_CRISPRroots-1.3>"
```

We recommend to first run Snakemake with `--dryrun`.  
This displays what will be done without executing it and highlights if any input file is missing. 
In the commands above, `--cores` specifies the maximum number of cores used in parallel by Snakemake.
The option `--singularity-args` is necessary to provide snakemake all the parameters to be passed to the singularity run. These are usually the `--bind` parameters that link local folders to the singularity.

The pre-computed results/reports for the test dataset are available in the folder 
*CRISPRroots_test_dataset/QPRT_DEL268T_chr16_10M-40M/pre-computed*.

Hint: Snakemake allows to visualize the jobs as a graph (directed acyclic graph, or DAG), 
highlighting the jobs completed and those to be run in different ways. 
To create an *svg* plot of your DAG, run the following command:
```shell
    snakemake -s <path_to_CRISPRroots>/run.smk --dag | dot -Tsvg > dag.svg
```
To learn more about the DAG and the visualization of jobs please visit the Snakemake tutorial
at https://snakemake.readthedocs.io/en/stable/tutorial/basics.html (*Step 4: Indexing read alignments and visualizing the DAG of jobs*).

### Advanced usage

The pipeline can also be used to accomplish only specific tasks. 
For example, to only perform the the pre-processing of the dataset and read mapping, 
you can run CRISPRroots with the rule flag `preproc_and_map` at the end:
```shell
$ snakemake -s <path_to_CRISPRroots>/run.smk --cores <int> --use-singularity --singularity-args
"--bind <global_path_to_QPRT_DEL268T_chr16_10M-40M>:<global_path_to_QPRT_DEL268T_chr16_10M-40M>
--bind <global_path_to_resources>:<global_path_to_resources>
--bind <global_path_to_CRISPRroots-1.3>:<global_path_to_CRISPRroots-1.3>" preproc_and_map
```
The rule preproc_and_map only executes the pre-processing and mapping steps. 

Relevant predefined target rules include:
*   **preproc_and_map**: Executes all pre-processing and mapping steps. Output:
    - Mapped reads: `<path_to_results_folder>/2_sortaligned/<sample name>/Aligned.Sorted.bam` 
    - Preprocessing and fastQC/multiQC quality reports: `<path_to_results_folder>/preproc/`
*   **variants_to_genome**: executes the rules necessary to produce files containing 
    filtered variants between each sample and the reference genome. 
    Output: `<path_to_results_folder>/6_GATK_variants/<sample name>/variants_filtered.vcf`
*   **on_target_check**: executes the on-target editing assessment. 
    Output: `<path_to_report_folder>/on_target_knockin.xlsx`; 
    `<path_to_report_folder>/on_target_knockout.xlsx`
*   **get_variated_genome**: produces a variant-aware version of the reference genome, 
    in which variants discovered from the RNA-seq are introduced in the reference sequence. 
    Output: `<path_to_results_folder>/6_GATK_variants/variated_genome.fa`
*   **get_lib_type**: assesses the library type with RSeQC. 
    Output: `<path_to_results_folder>/2-1_RSeQC_libtype/`
*   **preproc_and_map**: runs the reads pre-processing and mapping.  
    Mapping output: `<path_to_results_folder>/2_sortaligned/`      
    Mapping statistics: `<path_to_report_folder>/report/mapping_stats.xlsx`  
    Pre-processing results: `<path_to_results_folder>/preproc/`  
    Pre-processing statistics: `<path_to_report_folder>/report/multiqc_samples_stats.xlsx` 

NB: In the test dataset, `<path_to_results_folder>` and `<path_to_report_folder>` correspond to the 
subfolders *results* and *report* in that will be generated inside 
*CRISPRroots_test_dataset/QPRT_DEL268T_chr16_10M-40M/* after completing the pipeline.

A useful flag for the execution of specific tasks is `--notemp`. 
This avoids removing output files defined as temporary in the pipeline 
(e.g. partially processed reads). It can be convenient to use it when executing only a 
part of the pipeline, to avoid the removal of temporary files that will need to be recreated 
if required by a subsequent execution of the pipeline.

### Running the pipeline in a computer cluster

An example of how to set up `CRISPRroots` to run it with the Slurm Workload Manager is given in the test directory as `cluster_run.sh`. Configurations are given in `cluster_config.yaml`. To test the `cluster_run.sh` script on your system, first you need to fill in the fields:
* `<set_global_path_to_run.smk>` : path to the run.smk script in the CRISPRroots-1.3 folder
* `--singularity-args "<set_singularity_args_eg._binds>"` : all your singularity arguments,
such as bindings to the folders containing data, resources, and code

Once the script is complete, you can run it simply as:
```shell
  ./cluster_run.sh
```
You can add the name of a target rule to run only a part of the pipeline
```shell
  ./cluster_run.sh [target_rule]
```
      
## Output files
The pipeline’s output files are collected in two folders: **report** and **results**. 

The **report** folder contains the main output, including the candidate off-targets and the knockin/knockout
assessment. Results regarding differential expression and processing statistics (data quality and mapping) are also present.

The **results** folder contains all the output files of the pipeline that are not marked as `temp` (use `--notemp` to persist all files).

Please refer to the *CRISPRroots_Manual.pdf* included in the CRISPRroots software folder for a complete description fo the 
output and its content.
## Copyright

Copyright 2021 by the contributors:

Giulia Corsi <giulia@rth.dk>, Veerendra Gadekar <veer@rth.dk>

GNU GENERAL PUBLIC LICENSE

This is a free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License, either version 3 of the License, or (at your option) any later version. 
See [http://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citations

If you use CRISPRroots in your publication please cite:  
**CRISPRroots: on- and off-target assessment of RNA-seq data in CRISPR-Cas9 edited cells**  
Corsi GI, Gadekar VP, Gorodkin J, Seemann SE. *Nucleic Acids Research (2021, in press)*.

## Contact

In case of problems or bug reports, please contact <software+crisprroots@rth.dk>
