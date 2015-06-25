import os
import sys
from reaktoro.core import *

###############################################################################
# The following is needed to ensure that PyYAML uses OrderedDict instead of
# regular dict. This is needed to preserve the order of the YAML elements.
# This workaround is given at:
# http://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-
# mappings-as-ordereddicts/21048064#21048064
###############################################################################
import yaml
from yaml.representer import Representer
from yaml.constructor import Constructor, MappingNode, ConstructorError
import collections

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)
###############################################################################

# The Database instance used by the ChemicalEditor instance
database = None

# The ChemicalEditor instance used to create ChemicalSystem and ReactionSystem instances
editor = None

# The ChemicalSystem instance
system = None

# The ReactionSystem instance
reactions = None

# The dictionary of ChemicalState instances
states = {}

# The file used to output the results


def parseNumberWithUnits(word, default_units):
    word = str(word)
    words = word.split()
    return (float(eval(word)), default_units) if len(words) == 1 \
        else (float(eval(words[0])), words[1])


def splitKeywordIdentifier(key):
    ispace = key.find(' ')
    return (key, None) if ispace == -1 else \
        (key[:ispace], key[ispace+1:])


def childrenWithKeyword(node, keyword):
    children = {}
    for key, child in node.iteritems():
        kwd, identifier = splitKeywordIdentifier(key)
        if kwd == keyword:
            identifier = len(children) if identifier is None else identifier
            children[identifier] = child
    return children


def processChemicalSystem(node, identifier):
    # Update a ChemicalEditor instance with an AqueousPhase instance
    def addAqueousPhase(node):
        phase = node.get('AqueousPhase')
        if phase is not None:
            species = phase.get('Species')
            editor.addAqueousPhase(species)

    # Update a ChemicalEditor instance with a GaseousPhase instance
    def addGaseousPhase(node):
        phase = node.get('GaseousPhase')
        if phase is not None:
            species = phase.get('Species')
            editor.addGaseousPhase(species)

    # Update a ChemicalEditor instance with a MineralPhase instances
    def addMineralPhases(node):
        phases = node.get('MineralPhases', [])
        phases = phases.split()
        for phase in phases:
            editor.addMineralPhase(phase)

    # Initialize the Database instance
    global database
    database = node.get('Database', 'supcrt98.xml')
    database = Database(database)

    # Initialize the ChemicalEditor instance
    global editor
    editor = ChemicalEditor(database)

    # Process the aqueous, gaseous and mineral phases
    addAqueousPhase(node)
    addGaseousPhase(node)
    addMineralPhases(node)

    # Initialize the ChemicalSystem instance
    global system
    system = ChemicalSystem(editor)

    print system

def processPhreeqc(node, identifier):
    # Get the Database and Input entries in the Phreeqc block
    database_filename = node.get('Database')
    input_filename = node.get('Input')

    # Assert the Database entry has been provided
    assert database_filename is not None, \
        'Expecting the `Database` entry in the `Phreeqc` block.'

    # Assert the Input entry has been provided
    assert input_filename is not None, \
        'Expecting the `Input` entry in the `Phreeqc` block.'

    # Create the Phreeqc instance
    phreeeqc = Phreeqx(database_filename, input_filename)

    # Set the global chemical system instance
    global system
    system = ChemicalSystem(phreeeqc)

    # Create a chemical state instance to hold the chemical state of Phreeqc
    state = ChemicalState(phreeeqc)

    print '--------------------------------------------------------------------'
    print 'Printing the resulting chemical system from the `Phreeqc` block...'
    print '--------------------------------------------------------------------'
    print system

    print '--------------------------------------------------------------------'
    print 'Printing the resulting chemical state from the `Phreeqc` block...'
    print '--------------------------------------------------------------------'
    print 'This chemical state can be referenced as `PhreeqcState`.'
    print '--------------------------------------------------------------------'
    print state

    # Check if a post-equilibration step was requested
    if node.get('Equilibrate', False) in ['True', 'true']:
        print '--------------------------------------------------------------------'
        print 'Performing the post-equilibration step as requested...'
        print '--------------------------------------------------------------------'
        equilibrate(state)
        print '--------------------------------------------------------------------'
        print 'Printing the new state of `PhreeqcState` after equilibration...'
        print '--------------------------------------------------------------------'
        print state

    # Store the Phreeqc state under the name `PhreeqcState`
    global states
    states['PhreeqcState'] = state


def processReactionSystem(value, identifier):

    # Update the ChemicalEditor instance with a MineralReaction instance
    def addMineralReaction(name, node):
        # Create a MineralReaction instance and configure it
        reaction = editor.addMineralReaction(name)

        # Get the equation of the reaction
        equation = node.get('Equation')

        # Assert the equation of the reaction was provided
        assert equation is not None, 'Expecting the `Equation` keyword in ' \
            'the `MineralReaction %s` block.' % name

        # Set the equation of the MineralReaction instance
        reaction.setEquation(equation)

        # Get the specific surface area of the mineral
        ssa = node.get('SpecificSurfaceArea')

        # Assert the specific surface area of the mineral was provided
        assert ssa is not None, 'Expecting the `SpecificSurfaceArea` ' \
            'keyword in the `MineralReaction %s` block.' % name

        # Set the specific surface area of the mineral
        ssa, units = parseNumberWithUnits(ssa, 'm2/g')
        reaction.setSpecificSurfaceArea(ssa, units)

        # Iterate over all blocks with keyword `Mechanism`
        mechanisms = childrenWithKeyword(node, 'Mechanism')
        for key, mechanism in mechanisms.iteritems():
            # Get the rate constant of the current mineral mechanism
            k0 = mechanism.get('RateConstant')

            # Assert the rate constant was provided
            assert k0 is not None, 'Expecting a `RateConstant` ' \
                'keyword in the `Mechanism %s` block inside of the ' \
                '`MineralReaction %s` block.' % (key, name)

            # Parse the rate constant for its value and units
            k0, k0_units = parseNumberWithUnits(k0, 'mol/(m2*s)')

            # Get the Arrhenius activation energy of the current mineral mechanism
            Ea = mechanism.get('ActivationEnergy')

            # Assert the activation energy was provided
            assert Ea is not None, 'Expecting a `ActivationEnergy` ' \
                'keyword in the `Mechanism %s` block inside of the ' \
                '`MineralReaction %s` block.' % (key, name)

            # Parse the activation energy for its value and units
            Ea, Ea_units = parseNumberWithUnits(Ea, 'J/mol')

            # Get the powers `p` and `q` for the mineral mechanism
            power_p = mechanism.get('PowerP')
            power_q = mechanism.get('PowerQ')

            # Collect all nodes with keywords `ActivityPower` and `PressurePower`
            activity_powers = childrenWithKeyword(mechanism, 'ActivityPower')
            pressure_powers = childrenWithKeyword(mechanism, 'PressurePower')

            # Create a list of MineralCatalyst instances for the current mineral mechanism
            catalysts = []
            for species, power in activity_powers.iteritems():
                catalysts.append(MineralCatalyst(species, 'activity', power))
            for species, power in pressure_powers.iteritems():
                catalysts.append(MineralCatalyst(species, 'pressure', power))

            # Create a MineralMechanism instance
            mechanism = MineralMechanism()
            mechanism.setRateConstant(k0, k0_units)
            mechanism.setActivationEnergy(Ea, Ea_units)
            if power_p is not None: mechanism.setPowerP(power_p)
            if power_q is not None: mechanism.setPowerQ(power_q)
            if catalysts != []: mechanism.setCatalysts(catalysts)

            # Finally add the current MineralMechanism to the MineralReaction instance
            reaction.addMechanism(mechanism)

    # Update a ChemicalEditor instance with MineralReaction instances
    def addMineralReactions(node):
        children = childrenWithKeyword(node, 'MineralReaction')
        for name, child in children.iteritems():
            addMineralReaction(name, child)

    # Assert the ChemicalSystem instance has been previously defined
    assert system is not None, 'Expecting a `ChemicalSystem` block before the ' \
        '`ReactionSystem` block.'

    # Process the mineral reactions
    addMineralReactions(value)

    # Initialize the ReactionSystem instance
    global reactions
    reactions = ReactionSystem(editor)


def processEquilibrium(node, identifier):

    # Assert the ChemicalSystem instance has been initialized before
    auxstr = ' ' + identifier if identifier != None else ''
    assert system is not None, 'A `ChemicalSystem` block must be ' \
        'defined before the `Equilibrium%s` block.' % auxstr

    # Assert the block Mixture was specified
    assert node.get('Mixture') != None, 'Expecting a Mixture block in the' \
        '`Equilibrium%s` block.' % auxstr

    # Create the EquilibriumProblem instance
    problem = EquilibriumProblem(system)

    # Create a default partition of the chemical system
    partition = Partition(system)

    # Initialize the ChemicalState instance
    state = ChemicalState(system)

    # The list of triplets (phase, volume, units) listed in ScaleVolume block
    scaled_phase_volumes = []

    # Process the Temperature block
    def processTemperature(node, identifier):
        temperature, units = parseNumberWithUnits(node, 'kelvin')
        problem.setTemperature(temperature, units)

    # Process the Pressure block
    def processPressure(node, identifier):
        pressure, units = parseNumberWithUnits(node, 'pascal')
        problem.setPressure(pressure, units)

    # Process the Mixture block
    def processMixture(node, identifier):
        for name, amount in node.iteritems():
            # Check if `name` points to some some chemical state in `states`
            if states.has_key(name):
                problem.add(states[name], amount)
            # Then `name` should be the name of a species or generic compound
            else:
                amount, units = parseNumberWithUnits(amount, 'mol')
                problem.add(name, amount, units)

    # Process the InertSpecies block
    def processInertSpecies(node, identifier):
        names = [name for name in node]
        partition.setInertSpecies(names)
        problem.setPartition(partition)
        for species, amount in node.iteritems():
            amount, units = parseNumberWithUnits(amount, 'mol')
            state.setSpeciesAmount(species, amount, units)

    # Process the ScaleVolume block
    def processScaleVolume(node, identifier):
        for phase, volume in node.iteritems():
            volume, units = parseNumberWithUnits(volume, 'm3')
            scaled_phase_volumes.append((phase, volume, units))

    print 'Processing `Equilibrium%s`...' % auxstr

    # Initialize the dictionary of functions that process the keyword blocks
    processors = {}
    processors['Temperature'] = processTemperature
    processors['Pressure'] = processPressure
    processors['Mixture'] = processMixture
    processors['InertSpecies'] = processInertSpecies
    processors['ScaleVolume'] = processScaleVolume

    # Iterate over all specified blocks
    for key, child in node.iteritems():
        keyword, child_id = splitKeywordIdentifier(key)
        processor = processors.get(keyword)
        assert processor is not None, 'Unknown keyword `%s` in ' \
            '`Equilibrium%s` block.' % (keyword, auxstr)
        processor(child, child_id)

    # Perform the equilibrium calculation
    res = equilibrate(state, problem)

    # Perform the scale of the phase volumes if any
    for phase, volume, units in scaled_phase_volumes:
        state.setPhaseVolume(phase, volume, units)

    # Store the calculate chemical state in a dictionary of chemical states
    states[identifier] = state

    print 'Successfully solved `Equilibrium%s` in %d iterations and %f seconds.' \
        % (auxstr, res.optimum.iterations, res.optimum.time)

    print state


def processPlots(plotnodes, plots):
    for plotnode, plot in zip(plotnodes.itervalues(), plots):
        x = plotnode.get('x')
        y = plotnode.get('y')
        xtitle = plotnode.get('xtitle')
        ytitle = plotnode.get('ytitle')
        legend = plotnode.get('legend')
        frequency = plotnode.get('frequency')
        if x is not None: plot.x(x)
        if y is not None: plot.y(y)
        if xtitle is not None: plot.xtitle(xtitle)
        if ytitle is not None: plot.xtitle(ytitle)
        if legend is not None: plot.legend(legend)
        if frequency is not None: plot.legend(frequency)


def processEquilibriumPath(value, identifier):
    # Get the identifiers of the `From` and `To` states
    state1_id = value.get('From')
    state2_id = value.get('To')

    # Check if both `From` and `To` states have been specified
    auxstr = ' ' + identifier if identifier != None else ''
    assert state1_id is not None, 'Expecting a `From` statement in the ' \
        '`EquilibriumPath%s` block.' % auxstr
    assert state2_id is not None, 'Expecting a `To` statement in the ' \
        '`EquilibriumPath%s` block.' % auxstr

    # Get the ChemicalState instances from the global `states`
    state1 = states.get(state1_id)
    state2 = states.get(state2_id)

    # Check if both `from` and `to` states have been calculate before
    assert state1 is not None, 'The chemical state with identifier %s ' \
        'has not been calculated yet. Ensure this is calculated before the ' \
        '`EquilibriumPath%s` block.' % (state1_id, auxstr)
    assert state2 is not None, 'The chemical state with identifier %s ' \
        'has not been calculated yet. Ensure this is calculated before the ' \
        '`EquilibriumPath%s` block.' % (state2_id, auxstr)

    # Assert the ChemicalSystem instance has been initialized before
    assert system is not None, 'A `ChemicalSystem` block must be ' \
        'defined before the `EquilibriumPath%s` block.' % auxstr

    # Initialize the EquilibriumPath instance
    path = EquilibriumPath(system)

    # Collect the nodes with `Plot` keyword
    plotnodes = childrenWithKeyword(value, 'Plot')

    # Create as many ChemicalPlot instances
    plots = path.plots(len(plotnodes))

    # Initialize the ChemicalPlot instances
    processPlots(plotnodes, plots)

    # Solve the equilibrium path problem
    path.solve(state1, state2)


def processKineticPath(value, identifier):
    # Assert the ReactionSystem instance has been initialized before
    auxstr = ' ' + identifier if identifier != None else ''
    assert reactions is not None, 'A `ReactionSystem` block must be ' \
        'defined before the `KineticPath%s` block.' % auxstr

    # Get the `From` and `To` times
    t1 = value.get('From')
    t2 = value.get('To')

    # Check if `From` and `To` statements have been provided
    auxstr = ' ' + identifier if identifier != None else ''
    assert t1 is not None, 'Expecting a `From` statement in the ' \
        '`KineticPath%s` block.' % auxstr
    assert t2 is not None, 'Expecting a `To` statement in the ' \
        '`KineticPath%s` block.' % auxstr

    # Parse the time variables to get their values and units
    t1, units1 = parseNumberWithUnits(t1, 'seconds')
    t2, units2 = parseNumberWithUnits(t2, 'seconds')

    # Convert t1 and t2 to seconds
    t1 = convert(t1, units1, 'seconds')
    t2 = convert(t2, units2, 'seconds')

    # Get the identifier of the initial condition chemical state
    state_id = value.get('InitialCondition')

    # Check if the initial condition has been specified
    assert state_id is not None, 'Expecting a `InitialCondition` statement ' \
        'in the `KineticPath%s` block.' % auxstr

    # Get the ChemicalState instances from the global `states`
    state = states.get(state_id)

    # Check if the initial condition chemical state has been calculate before
    assert state is not None, 'The chemical state with identifier %s ' \
        'has not been calculated yet. Ensure this is calculated before the ' \
        '`KineticPath%s` block.' % (state_id, identifier)

    # Get the kinetic species
    kinetic_species = value.get('KineticSpecies', '')

    # Initialize the KineticPath instance
    path = KineticPath(reactions)

    # Set the names of the kinetic species
    path.setPartition('kinetic = %s' % kinetic_species)

    # Collect the nodes with `Plot` keyword
    plotnodes = childrenWithKeyword(value, 'Plot')

    # Create as many ChemicalPlot instances
    plots = path.plots(len(plotnodes))

    # Initialize the ChemicalPlot instances
    processPlots(plotnodes, plots)

    # Solve the kinetic path problem
    path.solve(state, t1, t2, 's')

    print state


def interpret(inputfile, outputfile=None):
    # Parse the YAML input file
    doc = yaml.load(inputfile)

    # Set the global outputfile file
    global output
    output = sys.stdout if outputfile == None else outputfile

    processors = {}
    processors['ChemicalSystem'] = processChemicalSystem
    processors['ReactionSystem'] = processReactionSystem
    processors['Equilibrium'] = processEquilibrium
    processors['EquilibriumPath'] = processEquilibriumPath
    processors['KineticPath'] = processKineticPath
    processors['Phreeqc'] = processPhreeqc

    for key, value in doc.iteritems():
        keyword, identifier = splitKeywordIdentifier(key)
        processors[keyword](value, identifier)


def main():
    # Get the command-line arguments (exclude the first, which is the name of this file)
    argv = sys.argv[1:]

    # Assert there exist at least one argument, the inputfile file
    assert len(argv) != 0, \
        'Expecting at least a inputfile file as argument ' \
        '(e.g., reaktoro inputfile.yaml)'

    # Assert there are no more than two arguments, the inputfile and outputfile files
    assert len(argv) <= 2, \
        'Expecting at most two arguments, an inputfile file and the outputfile file name ' \
        '(e.g., reaktoro inputfile.yaml outputfile.txt)'

    # Initialize the inputfile file
    inputfile = argv[0]
    inputfile = file(inputfile)

    # Initialize the outputfile file
    outputfile = os.path.splitext(inputfile)[0] + '.txt' \
        if len(argv) < 2 else argv[1]
    outputfile = file(outputfile, 'w')

    # Interpret the inputfile file
    interpret(inputfile, outputfile)


if __name__ == '__main__':
    main()
