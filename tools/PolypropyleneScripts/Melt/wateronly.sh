#!/bin/bash -l

#SBATCH -t 100:00:00
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --ntasks-per-node=44
#SBATCH --cpus-per-task=1

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=brian.maloney1@ucdconnect.ie

# Specifies the jobname
#SBATCH --job-name=wateronly
#SBATCH --output=wateronly_run.out
#SBATCH --error=wateronly.err

module load gromacs/2024.5
source /opt/software/el9/gromacs/2024.5/bin/GMXRC

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!)
gmx_mpi grompp -f step8.1_nvt.mdp -c onlywater.gro -p wateronly.top -o wateronly.tpr
gmx_mpi mdrun -s wateronly.tpr -deffnm wateronly


