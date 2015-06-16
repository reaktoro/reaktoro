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


def processChemicalSystem(doc):
    # Get the ChemicalSystem block from the document
    node = doc.get('ChemicalSystem')

    # Assert a ChemicalSystem block exist in the document
    assert node is not None, 'Expecting a block ChemicalSystem where you \
        define your chemical system.'

    # Initialize the Database instance
    database = node.get('Database', 'supcrt98.xml')
    database = Database(database)

    # Initialize the ChemicalEditor instance
    editor = ChemicalEditor(database)

    # Process the aqueous, gaseous and mineral phases
    addAqueousPhase(node, editor)
    addGaseousPhase(node, editor)
    addMineralPhases(node, editor)

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


def processTemperature(parent):
    return processValueWithUnits(parent, 'Temperature', 298.15)


def processPressure(parent):
    return processValueWithUnits(parent, 'Pressure', 1.0e+5)


def processEquilibriumMix(parent, problem):
    mix = parent.get('Mix')
    assert mix is not None, 'A `Mix` block is expected for equilibrium calculations.'
    lines = mix.split('\n')
    for line in lines:
        words = line.split()
        if words is []:
            continue
        assert len(words) == 3, 'Expecting a line with (1) a compound name, \
            (2) an amount, and (3) the units of the amount. For example, H2O 1.0 kg'
        compound = words[0]
        amount = float(words[1])
        units = words[2]
        problem.add(compound, amount, units)


def processEquilibrium(node, identifier):
    temperature = processTemperature(node)
    pressure = processPressure(node)
    problem = EquilibriumProblem(system)
    problem.setTemperature(temperature)
    problem.setPressure(pressure)
    processEquilibriumMix(node, problem)
    state = ChemicalState(system)
    equilibrate(state, problem)
    states[identifier] = state


def processEquilibriumBlocks(doc):
    pass


def interpret(script):
    # Parse the YAML script string
    doc = yaml.load(script)

    # Initialize the ChemicalSystem instance
    processChemicalSystem(doc)


    print system


script = """
ChemicalSystem:
    Database: supcrt98.xml
    AqueousPhase:
        Species: H2O(l) H+ OH- Na+ Cl- Ca+2 Mg+2 CO2(aq) HCO3-
    MineralPhases: Calcite Dolomite Quartz

Equilibrium StateIC:
    Temperature: 100 celsius
    Pressure: 300 bar
    Mix: |
        H2O 1 kg
        NaCl 1 mol
        CaCO3 200 g
        MgCO3 2 g
        SiO2 10 g
"""

interpret(script)
