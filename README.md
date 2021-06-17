# CRISPRroots

CRISPRroots: **R**NA-seq based **o**n-target and **o**ff-**t**arget assessment of Cas9-mediated edit**s**

## Installation and Configuration

### Prerequisites

The CRISPRroots pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/).  
We recommend snakemake installation via conda so that `--use-conda` flag can be used when calling Snakemake to let it automatically handle all the dependencies of the pipeline.

*   conda : ≥4.8.5
*   snakemake : ≥5.3.0

The pipeline was tested in a x86 64 GNU/Linux environment with Ubuntu v.18.04.1 (or newer) installed. All other software requirements are satisfied by the Conda environments defined in Snakemake, except for the R packages.  
We did not find a suitable combination of R package versions available for installation via Conda, hence the required R packages are automatically installed via `0.0_install_R_pkgs.R` script in the `snakemake/scripts` folder.

### Setup

CRISPRroots (version 1.0) is available to download from [https://rth.dk/resources/crispr/](https://rth.dk/resources/crispr/).  
The test dataset can be downloaded separately from [https://rth.dk/resources/crispr/crisprroots/downloads/CRISPRroots_test_dataset.tar.gz](here).
After downloading the software and the test dataset and un-packing, we recommend to explore the directory structure of the test dataset.

*  CRISPRroots_test_dataset
    *   resources
    *   QPRT_DEL268T_chr16_10M-40M (sample directory)

In the sample directory there is the configuration file `config.yaml` where the path to the resources, samples and CRISPRroots snakefile has to be set. This file also contains all the parameters defined for the execution of various steps in the pipeline. A copy of this file is also there in the CRISPRroots software package, this file can be used as template to create the configuration file for your own dataset.

For a complete list of parameters, options, and usage examples for CRISPRroots please read the Manual included in the CRISPRroots software folder.

## Testing the software

Assuming you have installed the software in the prerequisites, and have downloaded the software and test dataset, un-packed and setup the `config.yaml` file in sample directory as discussed above, you can run the pipeline for the test dataset from the sample directory using the following command:

      snakemake -s [path_to_CRISPRroots]/run.smk --use-conda --dry-run
      snakemake -s [path_to_CRISPRroots]/run.smk --use-conda

We recommend to first run Snakemake with `--dryrun`. This displays what will be done without executing it and highlights if any input file is missing. Once the pipeline finishes running it will generate the results and reports inside the current directory. Just for example the pre-computed results and report are available in the folder called "pre-computed" inside the sample directory.

## Pipeline usage

The pipeline can also be used to accomplish only specific tasks. For example, to only perform the the pre-processing of the dataset and read mapping, you can run CRISPRroots with the rule flag `preproc_and_map` at the end:

      snakemake -s [path_to_CRISPRroots]/run.smk --use-conda preproc_and_map

Other main rules ready for use are:

*   **variants_to_genome**: executes all of the rules necessary to produce files containing filtered variants (vcf ) between each sample and the reference genome. Output in: [path_to_results_folder]/6_GATK_variants/[sample name]/variants_filtered.vcf
*   **eSNPKaryotyping**: Executes the R package eSNP-Karyotyping for the analysis of genome integrity from RNA-seq. The standard workflow is modified to employ reads mapped with STAR instead of TopHat2\. Output in: [path_to_report_folder]/eSNPKaryotyping/
*   **on_target_check**: executes the on-target editing assessment. Output in: [path_to_report_folder]/on_target_knockin.xlsx; [path_to_report_folder]/on_target_knockout.xlsx
*   **get_variated_genome**: produces a variant-aware version of the reference genome, in which variants discovered from the RNA-seq are introduced in the reference sequence. Output in: [path_to_results_folder]/6_GATK_variants/variated_genome.fa
*   **get_lib_type**: assesses the library type with RSeQC. Output in: [path_to_results_folder]/2-2_RSeQC_libtype/

A useful flag for the execution of CRISPRroots in Snakemake is `--notemp`. This avoids removing output files defined as temporary in the pipeline (e.g. partially processed reads). It can be convenient to use it when executing only a part of the pipeline, to avoid the removal of temporary files that will need to be recreated if required by a subsequent execution of the pipeline.

## Running the pipeline in a computer cluster

CRISPRroots can also be launched on a computer cluster. An example of how to set up CRISPRroots to run it with the Slurm Workload Manager is given in the test directory and can be used as:

      [path_to_CRISPRroots]/cluster_run.sh

Note that the path to `run.smk` is hard-coded in `cluster_run.sh` and might need to be modified for usage outside the test folder. You can add the name of a target rule to run only a part of the pipeline as below:

      [path_to_CRISPRroots]/cluster_run.smk [target_rule]

## Copyright

Copyright 2021 by the contributors (see AUTHORS file)

GNU GENERAL PUBLIC LICENSE

This is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, either version 3 of the License, or (at your option) any later version. You should have received a copy of the GNU General Public License along with RIsearch, see file COPYING. If not, see [http://www.gnu.org/licenses/](https://www.gnu.org/licenses/).

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citations

If you use CRISPRroots in your publication please cite:  
**CRISPRroots: on- and off-target assessment of RNA-seq data in CRISPR-Cas9 edited cells**  
Corsi GI, Gadekar VP, Gorodkin J, Seemann SE. _Submitted_.

## Contact

In case of problems or bug reports, please contact <software+crisprroots@rth.dk>
