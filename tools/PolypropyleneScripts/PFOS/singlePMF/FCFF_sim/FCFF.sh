#!/bin/bash -l

#SBATCH -o out-%j
#SBATCH -e err-%j

#SBATCH --job-name=FCFF_sim
#SBATCH -t 90:00:00
#SBATCH -N 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=20

# E-mail address of recipient
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brian.maloney1@ucdconnect.ie

module load gromacs/2024.5

gmx_mpi convert-tpr -s md200.tpr -nsteps 400000000 -o md200.tpr
gmx_mpi mdrun -deffnm md200 -cpi md200.cpt -nb gpu
