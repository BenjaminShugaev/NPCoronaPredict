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
#SBATCH --job-name=waterslab_sim
#SBATCH --output=waterslab_run.out
#SBATCH --error=waterslab_run.err

module load gromacs/2024.5
source /opt/software/el9/gromacs/2024.5/bin/GMXRC

# Change working directory to current directory
cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!) 
gmx_mpi grompp -f step6.1_minimization.mdp -o step6.1_minimization.tpr -c waterslab.gro -r waterslab.gro -p waterslab.top
gmx_mpi mdrun -v -deffnm step6.1_minimization
gmx_mpi grompp -f step6.2_npt.mdp -o waterslabnpt.tpr -c step6.1_minimization.gro -r waterslab.gro -n waterslab.ndx -p waterslab.top
gmx_mpi mdrun -v -deffnm waterslabnpt
gmx_mpi grompp -f step6.3_nvt.mdp -o waterslabnvt.tpr -c waterslabnpt.gro -r waterslab.gro -n waterslab.ndx -p waterslab.top
gmx_mpi mdrun -v -deffnm waterslabnvt


