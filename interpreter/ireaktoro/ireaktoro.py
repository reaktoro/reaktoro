import argparse, os, sys
from reaktoro import *
from tabulate import tabulate

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
output = None

# The minimum width of the output bars
minbarwidth = 100

def basedir():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(os.path.realpath(__file__))

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

    # Check if the provided database file exists, if not try a built-in one
#     if not os.path.isfile(database):
#         print 'The provided database `%s` could not be found in the current working directory.' %  database
#         print 'Checking if a built-in database exists with this name.'
#         exedir = basedir()
#         databasedir = os.path.join(exedir, 'databases')
#         databaselocal = os.path.join(databasedir, database)
#         if not os.path.isfile(databaselocal):
#             raise RuntimeError('The provided database `%s` does not exist.' % database)
#         print 'Successfully found a built-in database with the same name `%s`.'  % database
#         database = databaselocal

    # Finally create the Database instance
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

    print >>output, system

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
    phreeeqc = Phreeqc(database_filename, input_filename)

    # Set the global chemical system instance
    global system
    system = ChemicalSystem(phreeeqc)

    # Create a chemical state instance to hold the chemical state of Phreeqc
    state = ChemicalState(phreeeqc)

    print >>output, '--------------------------------------------------------------------'
    print >>output, 'Printing the resulting chemical system from the `Phreeqc` block...'
    print >>output, '--------------------------------------------------------------------'
    print >>output, system

    print >>output, '--------------------------------------------------------------------'
    print >>output, 'Printing the resulting chemical state from the `Phreeqc` block...'
    print >>output, '--------------------------------------------------------------------'
    print >>output, 'This chemical state can be referenced as `PhreeqcState`.'
    print >>output, '--------------------------------------------------------------------'
    print >>output, state

    # Check if a post-equilibration step was requested
    if node.get('Equilibrate', False) in ['True', 'true']:
        print >>output, '--------------------------------------------------------------------'
        print >>output, 'Performing the post-equilibration step as requested...'
        print >>output, '--------------------------------------------------------------------'
        equilibrate(state)
        print >>output, '--------------------------------------------------------------------'
        print >>output, 'Printing the new state of `PhreeqcState` after equilibration...'
        print >>output, '--------------------------------------------------------------------'
        print >>output, state

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
    state = EquilibriumState(system)

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
        for species, val in node.iteritems():
            val, units = parseNumberWithUnits(val, 'mol')
            if convertible(units, 'mol'):
                state.setSpeciesAmount(species, val, units)
            else:
                state.setSpeciesMass(species, val, units)
                

    # Process the ScaleVolume block
    def processScaleVolume(node, identifier):
        for phase, volume in node.iteritems():
            volume, units = parseNumberWithUnits(volume, 'm3')
            scaled_phase_volumes.append((phase, volume, units))

    print >>output, 'Processing `Equilibrium%s`...' % auxstr

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
        state.scalePhaseVolume(phase, volume, units)

    # Store the calculate chemical state in a dictionary of chemical states
    states[identifier] = state

    print >>output, 'Successfully solved `Equilibrium%s` in %d iterations and %f seconds.' \
        % (auxstr, res.optimum.iterations, res.optimum.time)

    print >>output, state


def processPlots(plotnodes, plots):
    for plotnode, plot in zip(plotnodes.itervalues(), plots):
        x = plotnode.get('x')
        y = plotnode.get('y')
        xlabel = plotnode.get('xlabel')
        ylabel = plotnode.get('ylabel')
        legend = plotnode.get('legend')
        legendposition = plotnode.get('legendposition')
        showlegend = plotnode.get('showlegend')
        frequency = plotnode.get('frequency')
        
        # Assert both `x` and `y` entries were provided
        assert x is not None, 'Expecting a `x` statement for the plot.'
        assert y is not None, 'Expecting a `y` statement for the plot.'
        
        # Ensure default values for xlabel and ylabel are properly set
        xlabel = x if xlabel is None else xlabel
        ylabel = y if ylabel is None and type(y) is str else ylabel
        
        legend = y if legend is None else legend
        
        # Ensure both y and legend are lists
        if type(y) is not list: y = [y] 
        if type(legend) is not list: legend = [legend]
        
        assert len(legend) == len(y), 'Expecting the same number ' \
            'of entries in `y` and `legend` for the plot.'
         
        plot.xlabel(xlabel)
        if ylabel is not None: plot.ylabel(ylabel)
        
        plot.x(x)
        for label, quantity in zip(legend, y):
            plot.y(label, quantity)

        if legendposition is not None: plot.legend(legendposition)
        if showlegend is not None: plot.showlegend(showlegend)
        if frequency is not None: plot.frequency(frequency)
        
        xtics = plotnode.get('xtics')
        ytics = plotnode.get('ytics')
        if xtics: plot.xtics(str(xtics))
        if ytics: plot.ytics(str(ytics))


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
    state = KineticState(states.get(state_id))

    # Check if the initial condition chemical state has been calculate before
    assert state is not None, 'The chemical state with identifier %s ' \
        'has not been calculated yet. Ensure this is calculated before the ' \
        '`KineticPath%s` block.' % (state_id, identifier)

    # Get the kinetic species
    kinetic_species = value.get('KineticSpecies', [])

    # Initialize the KineticPath instance
    path = KineticPath(reactions)

    # Initialize the partition of the system
    partition = Partition(system)
    partition.setKineticSpecies(kinetic_species)
    
    # Set the names of the kinetic species
    path.setPartition(partition)

    # Collect the nodes with `Plot` keyword
    plotnodes = childrenWithKeyword(value, 'Plot')

    # Create as many ChemicalPlot instances
    plots = path.plots(len(plotnodes))

    # Initialize the ChemicalPlot instances
    processPlots(plotnodes, plots)

    # Solve the kinetic path problem
    path.solve(state, t1, t2, 's')

    print >>output, state


def processThermoProperties(node, identifier):
    # Get the path to the database file
    database = node.get('Database')

    # Ensure the Database entry was provided
    assert database != None, 'Expecting a `Database` entry in the ' \
        '`ThermoProperties` block (i.e., Database: supcrt98.xml).'

    # Get the temperature and pressure units
    tunits = node.get('TemperatureUnits', 'K')
    punits = node.get('PressureUnits', 'Pa')

    # Get the temperature and pressure points
    temperatures = node.get('Temperatures', [298.15]) # default to 298.15 K
    pressures = node.get('Pressures', [1e5]) # default to 1e5 Pa

    # Initialize the Database instance
    database = Database(database)

    # Initialize the Thermo instance
    thermo = Thermo(database)

    # Calculate the temperatures and pressures in standard units of K and Pa, respectively
    std_temperatures = [convert(val,  tunits, 'K') for val in temperatures]
    std_pressures = [convert(val,  punits, 'Pa') for val in pressures]

    # Get the list of reactions and species for which thermodynamic properties are calculated
    reactions = node.get('Reactions', [])
    species = node.get('Species', [])

    # Get the list of reaction and species properties to calculate
    reaction_properties = node.get('ReactionProperties', ['log(k)'])
    species_properties = node.get('SpeciesProperties', ['g'])

    # Transform the reaction and species properties to lower case
    reaction_properties = [x.lower() for x in reaction_properties]
    species_properties = [x.lower() for x in species_properties]

    # Get output options for the tables
    table_format = node.get('TableFormat', 'plain')
    alignment = node.get('Alignment', 'right')

    # Define a high-level function that returns the function to calculate a reaction property
    def reactionPropertyFunc(prop):
        if prop in ['logk', 'log(k)']:
            return thermo.logEquilibriumConstant
        if prop in ['lnk', 'ln(k)']:
            return thermo.lnEquilibriumConstant
        raise RuntimeError('The property `%s` specified in `ReactionProperties` is not supported.' % prop)

    # Define a high-level function that returns the function to calculate a species property
    def speciesPropertyFunc(prop):
        if prop in ['g', 'gibbsenergy']:
            return thermo.standardPartialMolarGibbsEnergy
        if prop in ['a', 'helmholtzenergy']:
            return thermo.standardPartialMolarHelmholtzEnergy
        if prop in ['u', 'internalenergy']:
            return thermo.standardPartialMolarInternalEnergy
        if prop in ['h', 'enthalpy']:
            return thermo.standardPartialMolarEnthalpy
        if prop in ['s', 'entropy']:
            return thermo.standardPartialMolarEntropy
        if prop in ['v', 'volume']:
            return thermo.standardPartialMolarVolume
        if prop in ['cp', 'heatcapacityconstp']:
            return thermo.standardPartialMolarHeatCapacityConstP
        if prop in ['cv', 'heatcapacityconstv']:
            return thermo.standardPartialMolarHeatCapacityConstV
        raise RuntimeError('The property `%s` specified in `SpeciesProperties` is not supported.' % prop)

    # Define a function that returns a table of reaction properties
    def reactionPropertyTable(reaction, prop):
        f = reactionPropertyFunc(prop)
        return [[f(T, P, reaction).val for T in std_temperatures] for P in std_pressures]

    # Define a function that returns a table of species properties
    def speciesPropertyTable(species, prop):
        f = speciesPropertyFunc(prop)
        return [[f(T, P, species).val for T in std_temperatures] for P in std_pressures]

    # Define a dictionary of property names
    propdict = {}
    propdict['logk'] = propdict['log(k)']             = 'log(K)'
    propdict['lnk']  = propdict['ln(k)']              = 'ln(K)'
    propdict['g']    = propdict['gibbsenergy']        = 'standard molar Gibbs energy (J/mol)'
    propdict['a']    = propdict['helmholtzenergy']    = 'standard molar Helmholtz energy (J/mol)'
    propdict['u']    = propdict['internalenergy']     = 'standard molar internal energy (J/mol)'
    propdict['h']    = propdict['enthalpy']           = 'standard molar enthalpy (J/mol)'
    propdict['s']    = propdict['entropy']            = 'standard molar entropy (J/(mol*K))'
    propdict['v']    = propdict['volume']             = 'standard molar volume (m3/mol)'
    propdict['cp']   = propdict['heatcapacityconstp'] = 'standard molar isobaric heat capacity (J/(mol*K))'
    propdict['cv']   = propdict['heatcapacityconstv'] = 'standard molar isochoric heat capacity (J/(mol*K))'

    # Define a function that outputs a formatted table of species/reaction properties
    def outputPropertyTable(table, entity, entitytype, prop):
        for (l, p) in zip(table, pressures): l.insert(0, p)
        headers = tuple(['P'] + ['T = %d' % x for x in temperatures])
        table = tabulate(table, headers=headers, tablefmt=table_format, numalign=alignment)
        caption = 'Calculated %s of %s %s' % (propdict[prop], entitytype, entity)
        width = max([len(x) for x in table.format().split('\n')])
        width = max(len(caption), width, minbarwidth)
        print >>output, '='*width
        print >>output, caption
        print >>output, '-'*width
        print >>output, table
        print >>output, '-'*width
        print >>output, 'TemperatureUnits: %s' % tunits
        print >>output, 'PressureUnits: %s' % punits
        print >>output, '='*width
        print >>output, ''

    # Output the calculated properties of the reactions
    for reaction in reactions:
        for prop in reaction_properties:
            table = reactionPropertyTable(reaction, prop)
            outputPropertyTable(table, reaction, 'reaction', prop)

    # Output the calculated properties of the species
    for sp in species:
        for prop in species_properties:
            table = speciesPropertyTable(sp, prop)
            outputPropertyTable(table, sp, 'species', prop)



# def processTransport(node, identifier):
#
#     # The Mesh instance
#     mesh = None
#
#     def processMesh(child, identifier):
#         dim = mesh.get('Dimension')
#
#         assert dim is not None, \
#             'The dimension of the mesh must be provided via keyword ' \
#             '`Dimension` in the `Mesh` block. ' \
#             'For example, Dimension: [1 km, 100 m].'
#
#         # Parse the components of the dimension vector in pairs (number, units)
#         dim = [parseNumberWithUnits(x, 'm') for x in dim]
#
#         # Convert the values to the same units
#         dim = [convert(x, units, 'm') for (x, units) in dim]
#
#         if len(dim) == 1:
#             mesh = ()
#
#         if len(dim) == 2:
#
#
#
#         discretization = mesh.get('Discretization')
#
#         assert discretization is not None, \
#             'The discretization of the mesh must be provided via keyword ' \
#             '`Discretization` in the `Mesh` block. ' \
#             'For example, Discretization: [100, 10].'
#
#         mesh = node.get('Mesh')
#
#     assert mesh is not None, \
#         'Expecting a `Mesh` block in the `Transport` block.'



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
    processors['ThermoProperties'] = processThermoProperties

    for key, value in doc.iteritems():
        keyword, identifier = splitKeywordIdentifier(key)
        processors[keyword](value, identifier)


def main():
    # Create a command-line argument parser
    parser = argparse.ArgumentParser(prog='Reaktoro')

    # Add the input argument
    parser.add_argument('input', type=str, \
        help='the relative path of the input file, including its name')

    # Add the output argument (optional)
    parser.add_argument('output', type=str, nargs='?', \
        help='the relative path of the output file, including its name')

    # Add the debug option (optional)
    parser.add_argument('-d', '--debug', action='store_true', \
        help='activates debug mode, which remotely communicates this ' \
            'application with the PyDev debugger')

    # Parse the command-line arguments (remove the first argument, which is the name of this file
    args = parser.parse_args(sys.argv[1:])

    # Check if debug is activated, and use pydev to start the remote debugger
    if args.debug is True:
        import pydevd; pydevd.settrace()

    # Provide a default output file name if none was provided
    if args.output is None:
        resultsdir = os.path.join(os.getcwd(), 'results')
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        inputname = os.path.split(args.input)[1]
        args.output = os.path.splitext(inputname)[0] + '.txt'
        args.output = os.path.join(resultsdir, args.output)

    # Initialize the input file
    inputfile = open(args.input, 'r')

    # Initialize the outputfile file
    outputfile = open(args.output, 'w')

    # Interpret the inputfile file
    interpret(inputfile, outputfile)


if __name__ == '__main__':
    main()
