#!/bin/bash
#SBATCH --partition=bgmp          ## Partition (like a queue in PBS)
#SBATCH --job-name=SNEK           ## Job Name
#SBATCH --output=OUT.txt          ## File in which to store job output
#SBATCH --error=ERR.txt           ## File in which to store job error messages
#SBATCH --time=4-00:00:00         ## Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                 ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1       ## Number of tasks to be launched per Node
#SBATCH --account=bgmp            ## Account used for job
#SBATCH --cpus-per-task=40        ## number of cpus (cores) per task
#SBATCH --mem=179G

## Activate environment
conda activate snakemake

## Run snakemake
/usr/bin/time -v snakemake \
--cores 40 \
--use-conda




