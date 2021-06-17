#!/bin/bash
#module load anaconda3/1
mkdir -p cluster_logs
# this script must be called from the sample data directory. In case of the CRISPRroots_test_dataset, it is QPRT_DEL268T_chr16_10M-40M directory.
# Below the path must be set to the location of snakemake directory 
snakemake -s ../../CRISPRroots/run.smk \
          -j 50 \
          --rerun-incomplete \
          --use-conda \
          --latency-wait 60 \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch --job-name={cluster.job-name} -t {cluster.time} -p {cluster.partition} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --output={cluster.output}" $1



