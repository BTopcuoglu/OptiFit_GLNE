#!/bin/bash

#SBATCH --job-name=smk_optifit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6GB
#SBATCH --time=5-00:00:00
#SBATCH --account=pschloss1
#SBATCH --partition=standard
#SBATCH --mail-user=armourc@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=%x-%j.out
#SBATCH --export=ALL

# List compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
	echo "Running on"
	scontrol show hostnames $SLURM_JOB_NODELIST
	echo -e "\n"
fi

# Connecting node to internet
source /etc/profile.d/http_proxy.sh

#activate snakemake
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glne

# Making output dir for snakemake cluster logs
mkdir -p logs/slurm/

# Initiating snakemake and running workflow in cluster mode
time snakemake --profile config/slurm/ --latency-wait 90 --rerun-incomplete 
