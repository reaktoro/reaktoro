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


def processValueWithUnits(parent, name, default):
    value = parent.get(name, default)
    if type(value) is str:
        words = value.split()
        value = float(words[0])
        if len(words) > 1:
            units = words[1]
            value = convert(value, units, 'kelvin')
    return value


def splitNumberUnits(word):
    assert type(word) in [int, float, str], \
        'Expecting a number with or without units.'
    if type(word) is str:
        ispace = word.find(' ')
        return (float(word), None) if ispace == -1 else \
            (float(word[:ispace]), word[ispace+1:])
    else:
        return (word, None)


def splitKeywordIdentifier(key):
    ispace = key.find(' ')
    return (key, None) if ispace == -1 else \
        (key[:ispace], key[ispace+1:])


def getTemperature(parent):
    value = parent.get('Temperature', 298.15)
    value, units = splitNumberUnits(value)
    return value if units is None else \
        convert(value, units, 'kelvin')


def getPressure(parent):
    value = parent.get('Pressure', 1.0e+5)
    value, units = splitNumberUnits(value)
    return value if units is None else \
        convert(value, units, 'pascal')


def getMixture(parent):
    mix = parent.get('Mix')
    mixture = []
    lines = mix.split('\n')
    for line in lines:
        words = line.split()
        if words == []:
            continue
        assert len(words) == 3, \
            'Expecting a line with (1) a compound name, ' \
            '(2) an amount, and (3) the units of the amount. ' \
            'For example, `H2O 1.0 kg`'
        compound = words[0]
        amount = float(words[1])
        units = words[2]
        mixture.append((compound, amount, units))
    return mixture


def processEquilibrium(value, identifier):
    # Assert the ChemicalSystem instance is not None
    assert system is not None, \
        'A ChemicalSystem block must be defined \
            before an Equilibrium block.'
    temperature = getTemperature(value)
    pressure = getPressure(value)
    mixture = getMixture(value)
    assert mixture != [], \
        'Expecting a `Mix` block in the Equilibrium block.'
    problem = EquilibriumProblem(system)
    problem.setTemperature(temperature)
    problem.setPressure(pressure)
    for compound, amount, units in mixture:
        problem.add(compound, amount, units)
    state = ChemicalState(system)
    equilibrate(state, problem)
    states[identifier] = state


def interpret(script):
    # Parse the YAML script string
    doc = yaml.load(script)

    processors = {}
    processors['ChemicalSystem'] = processChemicalSystem
    processors['Equilibrium'] = processEquilibrium

    for key, value in doc.iteritems():
        keyword, identifier = splitKeywordIdentifier(key)
        processors[keyword](value, identifier)

