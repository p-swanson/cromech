import openmm as mm
from   openmm import app
from   openmm.unit import *
from openmm.unit.constants import * 
import parmed
import numpy as np
import pandas as pd
import random
from reducedstatedatareporter import ReducedStateDataReporter
############################
#This file contains the code for uh
# for making a system containig a user specified number of LJ chains
# each with a specified length (number of beads)
# Parts
# 1) Setup the parameters (potentials, temp. etc.)
# 2) make required files (pdb, psf, etc.)
# 3) equilibration with CPU
# 4) production with GPU
# 5) nice :)
############################
# from
#https://doi.org/10.3390/ijms23169322
name = 'pentamer_eq'
# reduced units go here
r_eps = 1*dimensionless #0.238 *kilocalories_per_mole
r_sig = 1*dimensionless #3.4*angstrom
r_mass = 1 #39.9*amu

r_temperature = 1 #boltzmann_constant

# then we convert so openmm is happy i think this is dumb but i am going to do any ways
kB = BOLTZMANN_CONSTANT_kB*AVOGADRO_CONSTANT_NA
#kB = boltzmann_constant*AVOGADRO_CONSTANT_NA
sig = r_sig#*angstrom
eps = r_eps#*kilojoule/mole
mass = r_mass#*amu
temperature = np.true_divide(r_temperature*r_eps,kB)*kelvin
time_step = np.sqrt(np.true_divide(np.sqrt(r_mass)*np.square(sig),eps))*femtoseconds
gamma = 1000.0*0.5*time_step/picoseconds
r_c = 2.5*sig
K = 555.5*eps/sig**2

bond_length = 0.97*sig
box_size = 2000*sig
chain_len = 5
n_chains = 2000
natom    = chain_len*n_chains
#### iso config parameters
n_iso_run = 3
box_vecs = box_size*np.eye(3)
coords = []
################################ build system
from parmed import Structure, Atom, Bond, Angle, Dihedral

s  = Structure()
# 25 chains, one for each letter of the alphabet except Z so we have 5x5 grid
#chain_name = [f'{chr(i)}{j}' for j in range(5) for i in range(ord('A'),ord('Y')+1)]
chain_name = [f'{chr(i)}{j}' for j in range(80) for i in range(ord('A'),ord('Y')+1)]
resid_count = int()
total_count = int()
for n in chain_name:
  #resid_count = total_count
  for i in range(chain_len):
    s.add_atom(atom     =  Atom(name='X',  mass=mass),
              resname  = "LJ",
              resnum   =  resid_count,
              chain    = f'{n}')
    resid_count += 1
  #total_count += resid_count
for _ in range(n_chains):
      shift = random.randint(0,box_size)
      if i%2 == 0:
        plot_y = random.randint(0,box_size)
        plot_z = random.randint(0,box_size)
        for x in range(chain_len):
          coords.append((x+shift,plot_y,plot_z)) #for x in range(chain_len))
      else:
        plot_x1 = random.randint(0,box_size)
        plot_y1 = random.randint(0,box_size)
        for z in range(chain_len):
          coords.append((plot_x1,plot_y1,z+shift)) #for z in range(chain_len))

#s.positions = [(x, y, z) for x in range(int(n_chains**0.5)) for y in range(int(n_chains**0.5)) for z in range(chain_len)]
s.positions = coords
print('coords set')

bond_count = int()

for e in range(n_chains):
  #bond_count = total_bond
  for i in range(chain_len-1):
    s.bonds.append(Bond(s.atoms[i], s.atoms[i+1]))
    bond_count += 1
  #total_bond += bond_count
  bond_count += 1

s.save(f'{name}.pdb', overwrite=True)
s.save(f'{name}.psf', overwrite=True)
s.topology.setPeriodicBoxVectors(box_vecs)
########################################## add particles to the system
system = mm.System()

for i in range(natom):
  system.addParticle(s.atoms[i].mass)

system.setDefaultPeriodicBoxVectors(*box_vecs )
################################################
# epsilon = 0.8368 for all particles just save time by ading as a global parameter
energy_expression = '4*eps*((sig/r)^12-(sig/r)^6)'
force = mm.CustomNonbondedForce(energy_expression+';eps=sqrt(eps1*eps2); sig=(sig1+sig2)/2') #Lorentz-Berthelot combining rule

force.addPerParticleParameter('eps')
force.addPerParticleParameter('sig')

# Particles are assigned properties in the same order as they appear in the System object
for _ in range(natom):
    force.addParticle([eps,sig])

# Set force cutoff parameters
force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
force.setCutoffDistance(r_c)  

# Add exclusions between bonded residues

force.createExclusionsFromBonds([(i+e*chain_len,i+e*chain_len+1) for e in range(n_chains) for i in range(chain_len-1) ], 1) # exlcude bonded particles from LJ

### Create Harmonic force
force2 = mm.HarmonicBondForce()
### Add bonds
bond_count = int()
total_bond = int()
for o in range(n_chains):
  #bond_count = total_bond
  for _ in range(chain_len-1):
    force2.addBond(bond_count, bond_count+1, bond_length,K) # std LJ: 1.5*sig
    bond_count += 1
  #total_bond += bond_count
  bond_count += 1

### Add a force to remove Center of Mass motion to prevent drifting of molecule for trial runs
#force3 = mm.CMMotionRemover()

# Added forces to system
system.addForce(force)
system.addForce(force2)
#system.addForce(force3)
print('system set')

########################################## minimize and equilibrate with CPU
def CustomLangevinIntegrator(temperature=298.0*kelvin, collision_rate=91.0/picoseconds, timestep=1.0*femtoseconds):
    # Compute constants.
    kT = kB * temperature
    gamma = collision_rate

    dimensions =3
    # Create a new custom integrator.
    integrator = mm.CustomIntegrator(timestep)

    #
    # If dimensions == 2, set up a dummy variable to remove z-axial velocities
    #
    if dimensions == 2:
        integrator.addPerDofVariable("dumv", 1.0)
        integrator.setPerDofVariableByName("dumv", [mm.Vec3(x=1.0, y=1.0, z=0.0)])
    #
    # Integrator initialization.
    #
    integrator.addComputePerDof("sigma", "sqrt(kT/m)")
    integrator.addGlobalVariable("kT", kT) # thermal energy
    integrator.addGlobalVariable("T", temperature) # temperature
    integrator.addGlobalVariable("b", np.exp(-gamma*timestep)) # velocity mixing parameter
    integrator.addPerDofVariable("sigma", 0) 
    integrator.addPerDofVariable("x1", 0) # position before application of constraints

    #
    # Allow context updating here.
    #
    integrator.addUpdateContextState();

    # 
    # Velocity perturbation.
    #
    integrator.addComputePerDof("v", "sqrt(b)*v + sqrt(1-b)*sigma*gaussian")
    integrator.addConstrainVelocities();
    
    #
    # Metropolized symplectic step.
    #
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")
    if dimensions == 2: # To get a 2D system, make z-velocities zero when moving x
        integrator.addComputePerDof("x", "x + v*dumv*dt")
    else:
        integrator.addComputePerDof("x", "x + v*dt")
    integrator.addComputePerDof("x1", "x")
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m + (x-x1)/dt")

    #
    # Velocity randomization
    #
    integrator.addComputePerDof("v", "sqrt(b)*v + sqrt(1-b)*sigma*gaussian")
    if dimensions == 2: # Remove the resulting z-velocities to get the correct Kinetic Energy
        integrator.addComputePerDof("v", "v*dumv")

    return integrator

min_int  = CustomLangevinIntegrator(temperature, gamma, time_step)

#min_int = mm.LangevinIntegrator(temperature,gamma,time_step)
min_plat = mm.Platform.getPlatformByName('CPU')
pdb = app.PDBFile(f'{name}.pdb')
#pdb.topology.setPeriodicBoxVectors(box_vecs)
mini_sim = app.Simulation(pdb.topology,system,min_int,min_plat)
mini_sim.context.setPositions(s.positions)
#mini_sim.loadState('lb_prod_output2.xml')
mini_sim.minimizeEnergy()
print('minimized successfully')
#minimized_positions = mini_sim.context.getState(getPositions=True).getPositions()
mini_sim.context.setVelocitiesToTemperature(temperature)

# Reporters
mini_sim.reporters = []
mini_sim.reporters.append(app.XTCReporter(f'{name}_em.xtc', 500))
"""
mini_sim.reporters.append(app.StateDataReporter(f"{name}_em.csv", 500,
                                              time=True,
                                              density=True,
                                              potentialEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              volume=True))

"""
mini_sim.reporters.append(ReducedStateDataReporter(f"{name}_em.csv",500,3,r_temperature,step=True,potentialEnergy=True,kineticEnergy=True,temperature=True))
mini_sim.step(10000)
positions = mini_sim.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(mini_sim.topology, positions, open(f'{name}_em.pdb', 'w'))
mini_sim.saveState(f'{name}_em.xml') # save output config for equilbration run


##################################### production with GPU
"""                         
integrator = mm.CompoundIntegrator()
integrator.addIntegrator(mm.LangevinIntegrator(temperature, 1/picosecond, 2*femtoseconds))
integrator.addIntegrator(mm.VerletIntegrator(0.1*femtoseconds)) #

integrator = mm.LangevinIntegrator(temperature,1/picosecond,2*femtoseconds)

platform = mm.Platform.getPlatformByName('CUDA')
prop = {'CudaDeviceIndex': '0,1,2,3', 'CudaPrecision': 'single'}
simulation = app.Simulation(s.topology,
                            system,
                            integrator,platform,prop)

#simulation.context.setPositions(positions)
#simulation.minimizeEnergy()# minimize again just in case values too large for gpu
#integrator.setCurrentIntegrator(0)
# - Initialize velocities with random values at RT.
#simulation.loadState(f'{name}_em.xml')
simulation.loadState('mini_lb_to_sm_test_prod_output.xml')
simulation.context.setVelocitiesToTemperature(temperature)

# Reporters
simulation.reporters = []
simulation.reporters.append(app.DCDReporter(f'{name}_traj.dcd', 50000,enforcePeriodicBox=True))
#simulation.reporters.append(app.PDBReporter('ljtraj.pdb', 100),enforcePeriodicBox=True))
simulation.reporters.append(app.StateDataReporter(f"{name}_scalars.csv", 50000,
                                              time=True,
                                              density=True,
                                              potentialEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              volume=True))

simulation.step(5000000)
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(f'{name}_output.pdb', 'w'))
simulation.saveState(f'{name}_output.xml') # save output config for iso-config ensemble
"""
