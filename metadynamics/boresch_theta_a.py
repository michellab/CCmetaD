from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import mdtraj as md
import numpy as np

pdb_file = 'cc-di_fromtleap.pdb'
prmtop_file = 'cc-di.prmtop'
inpcrd_file = 'cc-di.inpcrd'

prot = md.load_pdb(pdb_file)

chain1 = [atom.index for atom in prot.topology.chain(0).atoms] #chain1
chain2 = [atom.index for atom in prot.topology.chain(1).atoms] #chain2
assembly = [atom.index for atom in prot.topology.atoms] #full assembly

coord = md.load(inpcrd_file, top = prmtop_file)

phi_array = md.compute_phi(coord)
psi_array = md.compute_psi(coord)

phi_atoms = phi_array[0].tolist()
psi_atoms = psi_array[0].tolist()

cas = [atom.index for atom in prot.topology.atoms if (atom.name == 'CA')]

residues= [str(res) for res in prot.topology.residues]

for res in residues:
    if res.startswith('ACE'):
        residues.remove(res)

for res in residues:
    if res.startswith('NHE'):
        residues.remove(res)

# anchor points for Boresch restraints
idx = 0

for i in residues:
    if i == 'ILE11':
        c_place = idx
    elif i == 'LYS16':
        b_place = idx
    elif i == 'ASN18':
        a_place = idx
    elif i == 'ASN52':
        A_place = idx
    elif i == 'TRP57':
        B_place = idx
    elif i == 'ILE59':
        C_place = idx
    idx+=1

c_idx = cas[c_place]
b_idx = cas[b_place]
a_idx = cas[a_place]
A_idx = cas[A_place]
B_idx = cas[B_place]
C_idx = cas[C_place]

# average alpha helix phi and psi dihedrals
phi_avg = (-57.0*np.pi)/180.0
psi_avg = (-47.0*np.pi)/180.0

# Ramachandran plot region for right-handed alpha-helix (in deg)
phi_min = -130.0
phi_max = -30.0
psi_min = -68.0
psi_max = 30.0

phi0 = (phi_min +phi_max)/2.0
psi0 = (psi_min+psi_max)/2.0
phicutoff = abs(phi_max-phi_min)/2.0
psicutoff = abs(psi_max-psi_min)/2.0

# defining the biasing range for the CV
cv_min = # in rad
cv_max = # in rad

# spring constants for the restraints
k_dihed = 100.0*kilojoules_per_mole


prmtop = AmberPrmtopFile(prmtop_file)
inpcrd = AmberInpcrdFile(inpcrd_file)
pdb = PDBFile(pdb_file)

system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds, hydrogenMass=1.5*amu, 
                             implicitSolvent=OBC2)

# definition of collective variable
boresch_theta_a = CustomAngleForce('theta')
boresch_theta_a.addAngle(b_idx, a_idx, A_idx)

bv = BiasVariable(boresch_theta_a, cv_min*radians, cv_max*radians, 0.1*radians)

# adding torsional restraints
phiforce = CustomTorsionForce('select({},{},{})'.format('step(-0.5*k_dihed*cos(theta-phi0)-(-0.5*k_dihed*cos_phicutoff))', '-0.5*k_dihed*cos(theta-phi0)', '-0.5*k_dihed*cos_phicutoff'))
phiforce.addGlobalParameter('k_dihed', k_dihed)
phiforce.addGlobalParameter('phi0', (phi0*np.pi/180.0)*radians)
phiforce.addGlobalParameter('cos_phicutoff', cos(phicutoff*np.pi/180.0))

psiforce = CustomTorsionForce('select({},{},{})'.format('step(-0.5*k_dihed*cos(theta-psi0)-(-0.5*k_dihed*cos_psicutoff))', '-0.5*k_dihed*cos(theta-psi0)', '-0.5*k_dihed*cos_psicutoff'))
psiforce.addGlobalParameter('k_dihed', k_dihed)
psiforce.addGlobalParameter('psi0', (psi0*np.pi/180.0)*radians)
psiforce.addGlobalParameter('cos_psicutoff', cos(psicutoff*np.pi/180.0))


for i in range(len(phi_atoms)):
    phiforce.addTorsion(phi_atoms[i][0], phi_atoms[i][1], phi_atoms[i][2], phi_atoms[i][3])
    psiforce.addTorsion(psi_atoms[i][0], psi_atoms[i][1], psi_atoms[i][2], psi_atoms[i][3])
    
system.addForce(phiforce)
system.addForce(psiforce)

for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)


meta = Metadynamics(system, [bv], 298.15*kelvin, 20.0, 0.8*kilojoule_per_mole, 100)

integrator = LangevinMiddleIntegrator(298.15*kelvin, 1/picosecond, 0.004*picoseconds)

platform = Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}

simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

save_freq = 1000

simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', save_freq))
simulation.reporters.append(StateDataReporter(stdout, save_freq, step=True,
        potentialEnergy=True, temperature=True))

steps = 1000000

meta.step(simulation, steps)

traj = md.load('output.dcd', top = prmtop_file)

cv_tuple=[]

for crd in traj.xyz: # getting atom coordinates in each frame to make the simulation context
    simulation.context.setPositions(crd)
    cv_tuple.append(meta.getCollectiveVariables(simulation)) # get CV per frame

fe = []
for i in meta.getFreeEnergy():
    fe.append(i._value)

fe_minimum = min(fe)
fe_norm = []
for i in range(len(fe)):
        fe_norm.append((fe[i]-fe_minimum)/4.184) # convert to kcal/mol as OpenMM uses kJ/mol
    
f = open('cv.txt', 'w')
for i in cv_tuple:
    f.write(str(i[0])+'\n') 
f.close()

f2 = open('fe_norm_kcal.txt', 'w')
for i in fe_norm:
    f2.write(str(i)+'\n')
f2.close()

# write out CV range
cv_range = np.arange(cv_min, cv_max, (cv_max-cv_min)/meta._widths[0])

with open('cv_range.txt', 'w') as cv_file:
    for cv_val in cv_range:
        cv_file.write(f"{cv_val}\n")


