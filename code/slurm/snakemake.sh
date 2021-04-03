#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=snakemake

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6GB
#SBATCH --time=5-00:00:00

# Account
#SBATCH --account=pschloss1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=armourc@umich.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=%x-%j.out

# Environment
#SBATCH --export=ALL

# Connecting node to internet
source /etc/profile.d/http_proxy.sh 

# List compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
	echo "Running on"
	scontrol show hostnames $SLURM_JOB_NODELIST
	echo -e "\n"
fi



#####################
#                   #
#  2) Job Commands  #
#                   #
#####################

# Making output dir for snakemake cluster logs
mkdir -p logs/slurm/

# Initiating snakemake and running workflow in cluster mode
snakemake --use-conda --verbose --profile config/slurm/  --latency-wait 90
