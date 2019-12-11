#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=mothurContigs.sh

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8000mb
#SBATCH --time=20:00:00

# Account
#SBATCH --account=pschloss
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=begumtop@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL

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
bash code/bash/getFiles.sh

bash code/bash/mothurReferences.sh

bash code/bash/mothurContigs.sh

bash code/bash/mothurLOO.sh

bash code/bash/mothurOptiFit.sh
