# Replicate CRISPRroots manuscript results
This folder contains the configuration files necessary to reproduce the results presened in the CRISPRroots manuscript. 
Each subfolder is dedicated to one of the 8 projects analyzed by CRISPRroots. 
The pipeline can be executed following the procedure described below. 

### Prerequisites

The CRISPRroots pipeline can be executed via 
[Snakemake](https://snakemake.readthedocs.io/en/stable/).  
To be executed in Sankemake, the pipeline requires:
* conda : ≥4.8.5
* snakemake : ≥5.32.00 

Additionally, you will need 
* python3

to run the automatic completion of the config files.

All other software requirements are satisfied by the Conda environments defined in Snakemake,
which are installed by starting Snakemake with the `--use-conda` flag.

The pipeline was tested in a x86 64 GNU/Linux environment with Ubuntu v.18.04.1 (or newer) 
installed.

### Resources and execution time
The pipeline was executed on a standard Linux node  with Intel&reg; Xeon&reg; CPU E5-2650, 60G RAM and 16 cores. 

The execution time varies largely between the different projects/datasets, depending mostly
on the number of samples and on the sequencing strategy or depth. The execution times
reported below were collected after running the full CRISPRroots pipeline, from raw sequencing data. 
The full pipeline includes the detection of germline variants from the RNA-seq data used for the definition of a new 
reference genome onto which potential off-targets with up to 6 mismatches are searched.

* from 6-8 hours for datasets with 3 replicates/condition (QPRT-INS/DEL, APOE, and PIK3CA-HET/HOM) 
* from 15-20 hours for datasets with 4 replicates/condition (GRIN2B and OGFOD1)

The variant calling (germline and somatic) were the most time-consuming steps, 
taking approx 5h each in the case of OGFOD1 and 11h (germline) and 5h (somatic) for GRIN2B-FW/REV.

Scheduling time of the queuing system is included in these estimates.
The time required by Conda to install all software requirements can vary from 10 minutes up to 1 hour. 
### RNA-seq data download
The RNA-seq data needs to be downloaded from https://www.ncbi.nlm.nih.gov/geo, 
from the following bioprojects: 
* QPRT-DEL and QPRT-INS: **GSE113734**
* GRIN2B_FW and GTIN2B_REV: **GSE114685**
* APOE: **GSE102956**
* PIK3CA-HOMO and PIK3CA-HET: **GSE126562**
* OGFOD1: **GSE130521**

The RNA-seq files in *fastq* format need to be saved under the respective folder, in the subfolder named ''samples'' and
named as described in the file ''samples_table.tsv'' of the corresponding project, with suffix *.fq.gz* for single end 
sequencing data and *_R1.fq.gz* or *_R2.fq.gz* for paired-end sequencing files. 
For instance, the correct sample names for the project QPRT-DEL are: QPRT_CTRL_1.fq.gz; QPRT_CTRL_2.fq.gz; QPRT_CTRL_3.fq.gz; QPRT_DEL268T_1.fq.gz; QPRT_DEL268T_2.fq.gz; QPRT_DEL268T_3.fq.gz.
### Completing the config file
The config file of each project comes pre-configured. You only need to run the script ''make_config.py'' included
in this folder to complete the config file by adding the global paths to your CRISPRroots, resources, and project folder.
From the current folder, run the following code after inserting your global path to the folders in ```<>```:
```shell
python3 make_config.py -c <path_to_CRISPRroots>/CRISPRroots/CRISPRroots-1.1 -d <path_to_CRISPRroots>/CRISPRroots/CRISPRroots-1.1/replicate_publication_results/QPRT-DEL -r <path_to_resources>/data/hg38_GRCh38_GENCODE33_FANTOMCAT
```
Repeat the same for all other project folders, substituting QPRT-DEL with each of
QPRT-INS, GRIN2B_FW, GTIN2B_REV, APOE, PIK3CA-HOMO, PIK3CA-HET, OGFOD1.

### Running the pipeline
You can run the pipeline by ```cd``` into one of the project folders and running CRISPRroots from there:
```shell
  cd QPRT-DEL
  snakemake -s <path_to_CRISPRroots>/run.smk --use-conda --dry-run 
  snakemake -s <path_to_CRISPRroots>/run.smk --use-conda --cores <int>
```

We recommend to first run Snakemake with `--dryrun`.  
This displays what will be done without executing it and highlights if any input file is missing. 
In the commands above, `--cores` specifies the maximum number of cores used in parallel by Snakemake.

Repeat the same procedure for each of the project folders, substituting QPRT-DEL with each of QPRT-INS, GRIN2B_FW, GTIN2B_REV, APOE, PIK3CA-HOMO, PIK3CA-HET, OGFOD1.

### Running the pipeline in a computer cluster

CRISPRroots can also be launched on a computer cluster. 
An example of how to set up a configuration file for CRISPRroots 
that allows it to run within the Slurm Workload Manager 
is provided in each of the project folders.
CRISPRroots can then be executed on the cluster as:
```shell
  cd QPRT-DEL
  ./<path_to_CRISPRroots>/cluster_run.sh
```
