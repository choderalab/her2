"""
Set up mTOR mutant simulations for SiegeTank.

@author John D. Chodera
@date 8 Aug 2014
"""

#
# IMPORTS
#

import os, os.path
import numpy

import pdbfixer

from simtk import openmm
from simtk import unit
from simtk.openmm import app

#
# PARAMETERS
#

# Path to put all output in
output_path = 'output'

# RCSB id
pdbid = '3PP0' # Pavletich structure of mTORdeltaN-mLST8 complex
chain_ids_to_keep = ['A'] # chains to retain
chain_id_to_mutate = 'A' # chain to mutate
pH = 7.0 # pH to model
keep_crystallographic_water = False # keep crystallographic waters?

# Single point mutants
point_mutants = ['L755S', 'V777L']

# DEBUG: Small problem for testing.
#pdbid = '1VII'
#chain_ids_to_keep = ['A']
#chain_id_to_mutate = 'F'
#point_mutants = ['F51A']

# Forcefield
ff_name = 'amber99sbildn'
water_name = 'tip3p'
padding = 11.0 * unit.angstroms
nonbonded_cutoff = 10.0 * unit.angstroms
nonbonded_method = app.CutoffPeriodic
max_minimization_iterations = 250
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
collision_rate = 5.0 / unit.picoseconds
barostat_frequency = 50
timestep = 1.0 * unit.femtoseconds
nsteps = 500 # number of steps to take for testing
nclones = 10 # number of clones per mutant

# Verbosity level
verbose = True

#===============================================================================
# DATA
#===============================================================================

three_letter_code = {
    'A' : 'ALA',
    'C' : 'CYS',
    'D' : 'ASP',
    'E' : 'GLU',
    'F' : 'PHE',
    'G' : 'GLY',
    'H' : 'HIS',
    'I' : 'ILE',
    'K' : 'LYS',
    'L' : 'LEU',
    'M' : 'MET',
    'N' : 'ASN',
    'P' : 'PRO',
    'Q' : 'GLN',
    'R' : 'ARG',
    'S' : 'SER',
    'T' : 'THR',
    'V' : 'VAL',
    'W' : 'TRP',
    'Y' : 'TYR'
}

one_letter_code = dict()
for one_letter in three_letter_code.keys():
    three_letter = three_letter_code[one_letter]
    one_letter_code[three_letter] = one_letter

def generate_pdbfixer_mutation_code(mutation):
    import re
    match = re.match('(\D)(\d+)(\D)', mutation)
    original_residue_name = three_letter_code[match.group(1)]
    residue_index = match.group(2)
    mutated_residue_name = three_letter_code[match.group(3)]
    return '%s-%s-%s' % (original_residue_name, residue_index, mutated_residue_name)

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

#
# Generate list of mutants.
#

npoint_mutants = len(point_mutants)

mutant_names = list()
mutant_codes = list()

# Append wild type (no mutation).
mutant_names.append('WT')
mutant_codes.append([])

# Append point mutants.
for i in range(npoint_mutants):
    mutation = point_mutants[i]
    mutant_names.append(mutation)
    mutant_codes.append([generate_pdbfixer_mutation_code(mutation)])

# Append all pairs of point mutants.
#for i in range(npoint_mutants):
#    for j in range(i+1, npoint_mutants):
#        mutation_i = point_mutants[i]
#        mutation_j = point_mutants[j]
#
#        mutant_names.append(mutation_i + '+' + mutation_j)
#        mutant_codes.append([generate_pdbfixer_mutation_code(mutation_i), generate_pdbfixer_mutation_code(mutation_j)])

#
# MAIN
#

# Create output directory.
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Open file to write all exceptions that occur during execution.
exception_filename = os.path.join(output_path, 'exceptions.out')
exception_outfile = open(exception_filename, 'w')

for (name, mutant) in zip(mutant_names, mutant_codes):
    print "%s : %s" % (name, str(mutant))

    try:
        # Initialize PDBFixer, retrieving PDB template.
        fixer = pdbfixer.PDBFixer(pdbid=pdbid)

        # Create directory to store files in.
        workdir = os.path.join(output_path, name)
        if not os.path.exists(workdir):
            os.makedirs(workdir)
            print "Creating path %s" % workdir

        # Hack to get chain id to chain number mapping.
        chain_id_list = [c.chain_id for c in fixer.structure.models[0].chains]

        # Build list of chains to remove
        chain_numbers_to_remove = list()
        for (chain_number, chain_id) in enumerate(chain_id_list):
            if chain_id not in chain_ids_to_keep:
                chain_numbers_to_remove.append(chain_number)

        # Remove all but desired chains.
        if verbose: print "Removing chains..."
        fixer.removeChains(chain_numbers_to_remove)

        if len(mutant) > 0:
            fixer.applyMutations(mutant, chain_id_to_mutate)

        if verbose: print "Adding missing atoms and residues..."
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH)
        fixer.removeHeterogens(keepWater=keep_crystallographic_water)

        # Write PDB file for solute only.
        if verbose: print "Writing pdbfixer output..."
        pdb_filename = os.path.join(workdir, 'pdbfixer.pdb')
        outfile = open(pdb_filename, 'w')
        app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
        outfile.close()

        # Refine in implicit solvent.
        if verbose: print "Minimizing energy in implicit solvent..."
        forcefield = app.ForceField(ff_name+'.xml', 'amber99_obc.xml')
        system = forcefield.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        if verbose: print "Creating simulation..."        
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        simulation = app.Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)
        simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
        state = simulation.context.getState(getEnergy=True)
        potential_energy = state.getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN after minimization.")
        if verbose: print "Potential energy after minimiziation: %.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)
        # Take a few steps to relax structure.
        if verbose: print "Taking a few steps..."
        simulation.step(nsteps)
        fixer.positions = simulation.context.getState(getPositions=True).getPositions()

        # Write PDB file for solute only.
        if verbose: print "Writing implicit simulation output..."
        pdb_filename = os.path.join(workdir, 'implicit.pdb')
        outfile = open(pdb_filename, 'w')
        app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
        outfile.close()

        # Solvate and create system.
        if verbose: print "Solvating..."
        forcefield = app.ForceField(ff_name+'.xml', water_name+'.xml')
        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(forcefield, padding=padding, model=water_name)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

        # Create simulation.
        if verbose: print "Creating simulation..."
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        # Write modeller positions.
        if verbose: print "Writing modeller output..."
        filename = os.path.join(workdir, 'modeller.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

        # Minimize energy.
        if verbose: print "Minimizing energy..."
        simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
        state = simulation.context.getState(getEnergy=True)
        potential_energy = state.getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN after minimization.")
        if verbose: print "Potential energy after minimiziation: %.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)

        # Write initial positions.
        filename = os.path.join(workdir, 'minimized.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

        # Assign temperature
        simulation.context.setVelocitiesToTemperature(temperature)

        # Take a few steps to relax structure.
        if verbose: print "Taking a few steps..."
        simulation.step(nsteps)

        # Write initial positions.
        if verbose: print "Writing positions..."
        filename = os.path.join(workdir, 'initial.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

        # Serialize to XML files.
        if verbose: print "Serializing to XML..."
        system_filename = os.path.join(workdir, 'system.xml')
        integrator_filename = os.path.join(workdir, 'integrator.xml')
        write_file(system_filename, openmm.XmlSerializer.serialize(system))
        write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))

        for clone_index in range(nclones):
            if verbose: print "Writing XML file for clone %d / %d..." % (clone_index+1, nclones)
            simulation.context.setVelocitiesToTemperature(temperature)
            state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
            state_filename = os.path.join(workdir, 'state%d.xml' % clone_index)
            serialized = openmm.XmlSerializer.serialize(state)
            write_file(state_filename, serialized)

    except Exception as e:
        print str(e)
        exception_outfile.write("%s : %s : %s\n" % (name, str(mutant), str(e)))
        exception_outfile.flush()

exception_outfile.close()

