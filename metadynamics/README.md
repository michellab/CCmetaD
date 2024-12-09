The **python scripts** provide examples of how to run sequential PMF calculation for the CC-Di coiled coil dimer with Boresch restraints. Appropriate ```cv_min``` and ```cv_max``` values must be chosen to bias each angle and torsion (in radians).

Note that each angle and torsion PMF calculation is done in the presence of only those restraints whose PMFs have been calculated beforehand (e.g. the theta_b PMF calculation is done with theta_a restrained to its equilibrium value, and the phi_a PMF calculation is done with both theta_a and theta_b restrained to their equilibrium values, etc.).

The notebook **metadynamics_analysis.ipynb** combines the simulation output files (cv_range.txt and fe_norm_kcal.txt) into a single file (requirement for free_energy_analysis.ipynb), and plots the PMF profiles. The files from each simulation replicate must be in directories called ```replicate-1/```, ```replicate-2/``` etc.

The notebook **free_energy_analysis.ipynb** analyses the free energy contributions from each PMF step. The individual steps and their replicates must be placed in directories called:
```
1_boresch_theta_a/
2_boresch_theta_b/
3_boresch_phi_a/
4_boresch_phi_b/
5_boresch_phi_c/
6_separation/
```

The directory **CC-Di_PMF_data" contains example data from the CC-Di simulations.

If RMSD restraints are used, their simulation files should be placed in directories called ```rmsd_bound/``` and ```rmsd_unbound/```.

An RMSD collective variable is defined as follows:
```
rmsd_cv = RMSDForce(pdb.positions, backbone_receptor)
```
where ```backbone_receptor``` is a list containing the atom indices of the receptor backbone atoms.

If RMSD restraints are used, the subsequent Boresch angle/torsion and separation PMF calculations must be done in the presence of RMSD restraints. A harmonic restraint can be placed on the RMSD CV as follows:
```
rmsd_cv = RMSDForce(pdb.positions, backbone_receptor) 
rmsd_receptor = CustomCVForce('0.5*k_rmsd*(rmsd-rmsd_target)^2')
rmsd_receptor.addCollectiveVariable('rmsd', rmsd_cv)
rmsd_receptor.addGlobalParameter('rmsd_target', 0.0)
rmsd_receptor.addGlobalParameter('k_rmsd', k_rmsd)
system.addForce(rmsd_receptor)
```
