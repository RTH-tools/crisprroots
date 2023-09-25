#!/bin/bash

mkdir -p cluster_logs
snakemake -s <set_global_path_to_run.smk> \
          -j 50 \
          --rerun-incomplete \
          --use-singularity \
          --latency-wait 60 \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch --job-name={cluster.job-name} -t {cluster.time} -p {cluster.partition} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --output={cluster.output}" \
	  --singularity-args "<set_singularity_args_eg._binds>" $1
