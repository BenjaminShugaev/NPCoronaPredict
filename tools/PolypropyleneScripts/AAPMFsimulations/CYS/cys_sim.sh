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
#SBATCH --job-name=cys_sim

module load gromacs/2024.5

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!) 
gmx_mpi convert-tpr -s md_50ns.tpr -nsteps 250000000 -o md_50ns.tpr
gmx_mpi mdrun -deffnm md_50ns -cpi md_50ns.cpt -nb gpu   
