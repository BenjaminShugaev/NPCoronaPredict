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
#SBATCH --job-name=ala_sim

module load gcc/10.1 gromacs/2023.3

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!) 
gmx_mpi convert-tpr -s md_50ns.tpr -nsteps 150000000 -o md_50ns.tpr
gmx_mpi mdrun -deffnm md_50ns -cpi md_50ns.cpt   
