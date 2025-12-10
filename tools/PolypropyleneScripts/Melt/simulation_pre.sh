gmx grompp -f step4.2_annealing.mdp -o ann.tpr -c eq.gro -r eq.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm ann
gmx grompp -f step4.3_nvt.mdp -o nvt.tpr -c ann.gro -r ann.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm nvt
gmx grompp -f step4.4_npt.mdp -o npt.tpr -c nvt.gro -r nvt.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm npt
