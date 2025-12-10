#!/bin/bash -l

#SBATCH -o out-%j
#SBATCH -e err-%j

#SBATCH -t 90:00:00
#SBATCH -N 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=20

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=brian.maloney1@ucdconnect.ie

# Specifies the jobname
#SBATCH --job-name=phe_sim


module load gromacs/2024.5

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!) 
gmx_mpi convert-tpr -s md_100ns.tpr -nsteps 250000000 -o md_100ns.tpr
gmx_mpi mdrun -deffnm md_100ns -cpi md_100ns.cpt -nb gpu   
