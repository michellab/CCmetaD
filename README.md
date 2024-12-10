# Coiled Coil metaDynamics

This repository contains data and code to run and analyse free energy calculations of coiled coils, according to the protocols in 'Assessment of the Topology and Oligomerisation States of Coiled
Coils Using Metadynamics with Conformational Restraints'.

## Contents

**designs**: scripts for design of coiled coils in ISAMBARD, coiled coil parameters from the design process and PDB structures used in simulations.

**metadynamics**: scripts for running and analysing metadynamics with conformational, orientational and positional restraints.

**restraints**: code for selecting anchor points and list of anchor points selected for each system along with equilibrium values.

## Requirements
[OpenMM](https://openmm.org/) v7.5

``conda install -c conda-forge openmm``

[MDTraj](https://mdtraj.org/1.9.4/index.html) v1.9.6

``conda install -c conda-forge mdtraj``
## Citation
TBD
