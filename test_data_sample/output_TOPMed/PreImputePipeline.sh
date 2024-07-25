#!/bin/bash

cd /DCEG/CGF/Bioinformatics/Production/Shilpa/Scripts/PreImputePipeline/test_data_sample/output_TOPMed
module load slurm
module load python3/3.10.2
DATE=$(date +%y%m%d)
mkdir -p logs_${DATE}
snakemake --cores=2 --unlock --configfile config.yaml
snakemake -pr --cluster "sbatch --partition=defq --cpus-per-task=10 --output=logs_${DATE}/snakejob_%j.out" --keep-going --rerun-incomplete --jobs 300 --latency-wait 120
