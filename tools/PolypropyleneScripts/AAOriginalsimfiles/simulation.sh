gmx grompp -f nvt.mdp -c system.gro -r system.gro -p system.top -n system.ndx -o nvt.tpr 
gmx mdrun -v -deffnm nvt 
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p system.top -n system.ndx -o npt.tpr 
gmx mdrun -v -deffnm npt 
gmx grompp -f pro.mdp -c npt.gro -t npt.cpt -p system.top -n system.ndx -o md_100ns.tpr 
gmx mdrun -v -deffnm md_100ns 
 
