A repo full of simulation files and notebooks describing how to simulate Polypropylene.

The **AAOriginalsim** folder contains files to create and run your own simulation boxes. The PMFgeneration notebook walks through the Gromacs commands and how everything is plotted. This folder also contains two simple scripts for PMF plotting, but these may need to be slightly altered to match your desired naming convention or file locations.

**AAPMFsimulations** contains the final frames and results involved with simulations we have previously ran. 

Inside the **Melt** folder, there are files used to run the polymer melt simulations. The orignal files were generated using CHARMM GUI polymer builder. The Melt notebook runs through the process of simulation.

The **PFOS** folder similarly runs through the simulations for PFOS, how to generate potentials and how to modify United Atom to examin no standard beads.

The UA notebook covers the basics of CoronaKMC simulations and scripts to plot results.

Finally, extra files such as the CalcLifschitzHamaker and UpscaleBeads are provided.These are provided with the standard UA, however also given here for completeness.
