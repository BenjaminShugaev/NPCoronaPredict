#!/bin/bash -l

# Set the number of nodes

#SBATCH -N 1

# Set the number of tasks/cores per node required 
#SBATCH -n 44

# Set the walltime of the job to 1 hour (format is hh:mm:ss)
#SBATCH -t 60:00:00

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=brian.maloney1@ucdconnect.ie

# Specifies the jobname
#SBATCH --job-name=ile_sim
#SBATCH --exclude=sonic60
#SBATCH --exclude=sonic61
#SBATCH --exclude=sonic56
#SBATCH --exclude=sonic74

module load gromacs/2024.5

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!) 
gmx_mpi convert-tpr -s md_100ns.tpr -nsteps 100000000 -o md_100ns.tpr
gmx_mpi mdrun -deffnm md_100ns -cpi md_100ns.cpt   
