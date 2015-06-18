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

# The ChemicalSystem instance
system = None

# The ReactionSystem instance
reactions = None

# The dictionary of ChemicalState instances
states = {}


def addAqueousPhase(doc, editor):
    phase = doc.get('AqueousPhase')
    if phase is not None:
        species = phase.get('Species')
        editor.addAqueousPhase(species)


def addGaseousPhase(doc, editor):
    phase = doc.get('GaseousPhase')
    if phase is not None:
        species = phase.get('Species')
        editor.addGaseousPhase(species)


def addMineralPhases(doc, editor):
    phases = doc.get('MineralPhases', [])
    phases = phases.split()
    for phase in phases:
        editor.addMineralPhase(phase)


def processChemicalSystem(value, identifier):
    # Initialize the Database instance
    database = value.get('Database', 'supcrt98.xml')
    database = Database(database)

    # Initialize the ChemicalEditor instance
    editor = ChemicalEditor(database)

    # Process the aqueous, gaseous and mineral phases
    addAqueousPhase(value, editor)
    addGaseousPhase(value, editor)
    addMineralPhases(value, editor)

    # Initialize the ChemicalSystem instance
    global system
    system = ChemicalSystem(editor)

    print system


def processValueWithUnits(parent, name, default):
    value = parent.get(name, default)
    if type(value) is str:
        words = value.split()
        value = float(words[0])
        if len(words) > 1:
            units = words[1]
            value = convert(value, units, 'kelvin')
    return value


def parseNumberWithUnits(word, default_units):
    word = str(word)
    words = word.split()
    return (float(word), default_units) if len(words) == 1 \
        else (float(words[0]), words[1])


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


def processEquilibrium(value, identifier):

    # Set the temperature of an EquilibriumProblem instance
    def setTemperature(problem):
        temperature = value.get('Temperature')
        if temperature != None:
            temperature, units = parseNumberWithUnits(temperature, 'kelvin')
            problem.setTemperature(temperature, units)

    # Set the pressure of an EquilibriumProblem instance
    def setPressure(problem):
        pressure = value.get('Pressure')
        if pressure != None:
            pressure, units = parseNumberWithUnits(pressure, 'pascal')
            problem.setPressure(pressure, units)

    # Set the mixture composition of an EquilibriumProblem instance
    def setMixture(problem):
        mixture = value.get('Mixture')
        assert mixture != None, 'Expecting a Mixture block in the' \
            '`Equilibrium %s` block.' % identifier
        for compound, amount in mixture.iteritems():
            amount, units = parseNumberWithUnits(amount, 'mol')
            problem.add(compound, amount, units)

    # Apply the scaling commands in the ScaleVolume block to the chemical state
    def applyScaleVolume(state):
        dic = value.get('ScaleVolume', {})
        for phase, volume in dic.iteritems():
            volume, units = parseNumberWithUnits(volume, 'm3')
            state.setPhaseVolume(phase, volume, units)

    print 'Processing Equilibrium %s...' % identifier

    # Assert the ChemicalSystem instance has been initialized before
    auxstr = ' ' + identifier if identifier != None else ''
    assert system is not None, 'A `ChemicalSystem` command block must be ' \
        'defined before the `Equilibrium%s` block.' % auxstr

    # Initialize the EquilibriumProblem instance
    problem = EquilibriumProblem(system)
    setTemperature(problem)
    setPressure(problem)
    setMixture(problem)

    # Initialize the ChemicalState instance
    state = ChemicalState(system)

    # Perform the equilibrium calculation
    res = equilibrate(state, problem)

    # Perform the scale of the phase volumes if any
    applyScaleVolume(state)

    # Store the calculate chemical state in a dictionary of chemical states
    states[identifier] = state

    print 'Successfully solved Equilibrium %s in %d iterations and %f seconds.' \
        % (identifier, res.optimum.iterations, res.optimum.time)

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
        '`EquilibriumPath%s` command block.' % auxstr
    assert state2_id is not None, 'Expecting a `To` statement in the ' \
        '`EquilibriumPath%s` command block.' % auxstr

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
    assert system is not None, 'A `ChemicalSystem` command block must be ' \
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
    # Get the `From` and `To` times
    t1 = value.get('From')
    t2 = value.get('To')

    # Check if `From` and `To` statements have been provided
    auxstr = ' ' + identifier if identifier != None else ''
    assert t1 is not None, 'Expecting a `From` statement in the ' \
        '`KineticPath%s` command block.' % auxstr
    assert t2 is not None, 'Expecting a `To` statement in the ' \
        '`KineticPath%s` command block.' % auxstr

    # Parse the time variables to get their values and units
    t1, units1 = parseNumberWithUnits(t1, 'seconds')
    t2, units2 = parseNumberWithUnits(t1, 'seconds')

    # Convert t1 and t2 to seconds
    t1 = convert(t1, units1, 'seconds')
    t2 = convert(t1, units2, 'seconds')

    # Get the identifier of the initial condition chemical state
    state_id = value.get('InitialCondition')

    # Check if the initial condition has been specified
    auxstr = ' ' + identifier if identifier != None else ''
    assert state_id is not None, 'Expecting a `InitialCondition` statement ' \
        'in the `KineticPath%s` command block.' % auxstr

    # Get the ChemicalState instances from the global `states`
    state = states.get(state_id)

    # Check if the initial condition chemical state has been calculate before
    assert state is not None, 'The chemical state with identifier %s ' \
        'has not been calculated yet. Ensure this is calculated before the ' \
        '`KineticPath%s` block.' % (state_id, identifier)

    # Assert the ReactionSystem instance has been initialized before
    assert reactions is not None, 'A `Reactions` command block must be ' \
        'defined before the `KineticPath%s` block.' % auxstr

    # Initialize the KineticPath instance
    path = KineticPath(reactions)

    # Collect the nodes with `Plot` keyword
    plotnodes = childrenWithKeyword(value, 'Plot')

    # Create as many ChemicalPlot instances
    plots = path.plots(len(plotnodes))

    # Initialize the ChemicalPlot instances
    processPlots(plotnodes, plots)

    # Solve the kinetic path problem
    path.solve(state, t1, t2)


def interpret(script):
    # Parse the YAML script string
    doc = yaml.load(script)

    processors = {}
    processors['ChemicalSystem'] = processChemicalSystem
    processors['Equilibrium'] = processEquilibrium
    processors['EquilibriumPath'] = processEquilibriumPath

    for key, value in doc.iteritems():
        keyword, identifier = splitKeywordIdentifier(key)
        processors[keyword](value, identifier)

