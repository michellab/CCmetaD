The script **separation.py** provides an example of how to run a separation PMF calculation for the CC-Di coiled coil dimer in the presence of Boresch restraints. The Boresch restraints are added to the simulation as follows:
```
boresch_theta_a = CustomAngleForce('0.5*k_boresch*(theta-theta_a_0)^2')
boresch_theta_a.addGlobalParameter('k_boresch', k_boresch)
boresch_theta_a.addGlobalParameter('theta_a_0', (theta_a_0*np.pi/180.0)*radians)
boresch_theta_a.addAngle(b_idx, a_idx, A_idx)
system.addForce(boresch_theta_a)

boresch_theta_b = CustomAngleForce('0.5*k_boresch*(theta-theta_b_0)^2')
boresch_theta_b.addGlobalParameter('k_boresch', k_boresch)
boresch_theta_b.addGlobalParameter('theta_b_0', (theta_b_0*np.pi/180.0)*radians)
boresch_theta_b.addAngle(a_idx, A_idx, B_idx)
system.addForce(boresch_theta_b)

boresch_phi_a = CustomTorsionForce('0.5*k_boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-phi_a_0); pi = 3.1415926535')
boresch_phi_a.addGlobalParameter('k_boresch', k_boresch)
boresch_phi_a.addGlobalParameter('phi_a_0', (phi_a_0*np.pi/180.0)*radians)
boresch_phi_a.addTorsion(c_idx, b_idx, a_idx, A_idx)
system.addForce(boresch_phi_a)

boresch_phi_b = CustomTorsionForce('0.5*k_boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-phi_b_0); pi = 3.1415926535')
boresch_phi_b.addGlobalParameter('k_boresch', k_boresch)
boresch_phi_b.addGlobalParameter('phi_b_0', (phi_b_0*np.pi/180.0)*radians)
boresch_phi_b.addTorsion(b_idx, a_idx, A_idx, B_idx)
system.addForce(boresch_phi_b)

boresch_phi_c = CustomTorsionForce('0.5*k_boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-phi_c_0); pi = 3.1415926535')
boresch_phi_c.addGlobalParameter('k_boresch', k_boresch)
boresch_phi_c.addGlobalParameter('phi_c_0', (phi_c_0*np.pi/180.0)*radians)
boresch_phi_c.addTorsion(a_idx, A_idx, B_idx, C_idx)
system.addForce(boresch_phi_c)
```

To run PMF calculations for each Boresch angle and torsion before the separation step, the CV in the script should be changed to one of the following depending on which angle/torsion is biased:
```
boresch_theta_a = CustomAngleForce('theta')
boresch_theta_a.addAngle(b_idx, a_idx, A_idx)
```
```
boresch_theta_b = CustomAngleForce('theta')
boresch_theta_b.addAngle(a_idx, A_idx, B_idx)
```
```
boresch_phi_a = CustomTorsionForce('theta')
boresch_phi_a.addTorsion(c_idx, b_idx, a_idx, A_idx)
```
```
boresch_phi_b = CustomTorsionForce('theta')
boresch_phi_b.addTorsion(b_idx, a_idx, A_idx, B_idx)
```
```
boresch_phi_c = CustomTorsionForce('theta')
boresch_phi_c.addTorsion(a_idx, A_idx, B_idx, C_idx)
```
Appropriate ```cv_min``` and ```cv_max``` values must be chosen to bias each angle and torsion (in radians).

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
If RMSD restraints are used, their simulation files should be placed in directories called ```rmsd_bound/``` and ```rmsd_unbound/```.

An RMSD collective variable is defined as follows:
```
rmsd_cv = RMSDForce(pdb.positions, backbone_receptor)
```
where ```backbone_receptor``` is a list containing the atom indices of the receptor backbone atoms.

If RMSD restraints are used, the subsequent Boresch angle and torsion PMF calculations must be done in the presence of RMSD restraints. A harmonic restraint can be placed on the RMSD CV as follows:
```
rmsd_cv = RMSDForce(pdb.positions, backbone_receptor) 
rmsd_receptor = CustomCVForce('0.5*k_rmsd*(rmsd-rmsd_target)^2')
rmsd_receptor.addCollectiveVariable('rmsd', rmsd_cv)
rmsd_receptor.addGlobalParameter('rmsd_target', 0.0)
rmsd_receptor.addGlobalParameter('k_rmsd', k_rmsd)
system.addForce(rmsd_receptor)
```
