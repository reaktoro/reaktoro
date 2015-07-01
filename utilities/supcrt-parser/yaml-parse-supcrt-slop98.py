#----------------------------------------------------------------------------------------------------
# Note:
# The following corrections are necessary in file slop07.dat
# before this script is used to parse it.
#
#  1) In species block AMORPHOUS-SILICA, correct the elemental formula to Si(1)O(2).
#  2) In species block LAURITE, change elemental formula to Ru(1)S(2).
#  3) In species block Ru(SO4)2-2, change elemental formula to Ru(1)S(2)O(8)-(2).
#  4) In species block Pd(OH)2(s), insert abbreviation Pd(OH)2.
#  5) In species block PdO(s), insert abbreviation PdO.
#  6) In species block H2O,g, correct the elemental formula from H(2)O to H(2)O(1).
#  7) In species blocks Pd+2, Rh+2, Rh+3, Ru+2, and Ru+3 remove the suffixes (II)ion or (III)ion
#     from their elemental formulas, which should be Pd(1), Rh(1), Rh(1), Ru(1) and Ru(1).
#  8) In species block Ce+4, correct the elemental formula to Ce(1) and formula to Ce.
#----------------------------------------------------------------------------------------------------
import yaml
from collections import OrderedDict
from complement import *

class SpeciesData:
    pass

# The list of names writen in capitals that should be converted to tilte style
capital_names = set(open('capital_names.txt').read().splitlines())

def molarMass(elemental_formula):
    return sum([atoms * elements[name] for name, atoms in elemental_formula])

def parseElementalFormula(elemental_formula):
    # Remove the trailing symbol (s) from the elemental formula of some minerals
    elemental_formula = elemental_formula.replace('(s)', '')

    # Split the elemental formula in words delimited by ( or )
    words = elemental_formula.replace(')', '(').split('(')
    words = [x for x in words if x != '']

    # Check if the elemental formula contains only one element
    if len(words) == 1:
        return [(words[0], 1)]

    # Collect the pairs (element, num_atoms) in the list elements
    elements = []
    for i in range(0, len(words), 2):
        if words[i] == '+' or words[i] == '-': # some formulas contain the suffix +(charge) or -(charge)
            break
        elements.append((words[i], int(words[i+1])))
    return elements

def jumpToSection(f, section):
    while(True):
        line = f.readline().strip()
        if line == section:
            f.readline() # skip one more line
            break

def parseGeneralData(f, data):
    # Parse the 1st line of the species data
    words = f.readline().strip().split()

    # Check if the first line is an empty line or a comment line
    if not words or words[0][0] == '*':
        return False

    data.name = words[0]
    data.formula = words[1]

    # Change the name of the species to title style if it is writen in caps
    if data.name in capital_names:
        data.name = data.name.title()

    # Parse the 2nd line of the species data
    words = f.readline().strip().split()

    data.abbreviation = words[0]
    data.elemental_formula = parseElementalFormula(words[1])

    # Parse the 3rd line of the species data
    words = f.readline().strip().split()
    data.references = words[0].replace('ref:', '').replace('REF:', '')
    data.date = words[1].title()
    return True

def parseAqueousData(f):
    # Create the species data (see Table 11 of SUPCRT92 paper)
    data = SpeciesData()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Parse the 4th line of the aqueous species data
    words = f.readline().strip().split()
    data.Gf = float(words[0])
    data.Hf = float(words[1])
    data.Sr = float(words[2])

    # Parse the 5th line of the aqueous species data
    words = f.readline().strip().split()
    data.a1 = float(words[0]) * 1.0e-01
    data.a2 = float(words[1]) * 1.0e+02
    data.a3 = float(words[2]) * 1.0e+00
    data.a4 = float(words[3]) * 1.0e+04

    # Parse the 6th line of the aqueous species data
    words = f.readline().strip().split()
    data.c1     = float(words[0]) * 1.0e+00
    data.c2     = float(words[1]) * 1.0e+04
    data.wref   = float(words[2]) * 1.0e+05
    data.charge = int(float(words[3]))

    # Set the type of the species
    data.type = 'Aqueous'

    return data

def parseGaseousData(f):
    # Create the species data (see Table 9 of SUPCRT92 paper)
    data = SpeciesData()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Parse the 4th line of the gaseous species data
    words = f.readline().strip().split()
    data.Gf = float(words[0])
    data.Hf = float(words[1])
    data.Sr = float(words[2])

    # Parse the 5th line of the gaseous species data
    words = f.readline().strip().split()
    data.a = float(words[0]) * 1.0e+00
    data.b = float(words[1]) * 1.0e-03
    data.c = float(words[2]) * 1.0e+05

    # Parse the 6th line of the gaseous species data
    words = f.readline().strip().split()
    data.Tmax = float(words[0])

    # Set the type of the species
    data.type = 'Gaseous'

    return data

def parseMineralData(f, nptrans):
    # Create the species data (see Table 10 of SUPCRT92 paper)
    data = SpeciesData()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Parse the 4th line of the mineral species data
    words = f.readline().strip().split()
    data.Gf = float(words[0])
    data.Hf = float(words[1])
    data.Sr = float(words[2])
    data.Vr = float(words[3])
    data.nptrans = nptrans

    # Parse the 5th line of the mineral species data
    if nptrans == 0:
        words = f.readline().strip().split()
        data.a = float(words[0]) * 1.0e+00
        data.b = float(words[1]) * 1.0e-03
        data.c = float(words[2]) * 1.0e+05
    else:
        data.a, data.b, data.c = [], [], []
        data.Ttr, data.Htr, data.Vtr, data.dPdTtr = [], [], [], []
        for i in range(nptrans+1):
            words = f.readline().strip().split()
            data.a.append(float(words[0]) * 1.0e+00)
            data.b.append(float(words[1]) * 1.0e-03)
            data.c.append(float(words[2]) * 1.0e+05)

            if i < nptrans:
                data.Ttr.append(float(words[3]))
                data.Htr.append(float(words[4]))
                data.Vtr.append(float(words[5]))
                data.dPdTtr.append(float(words[6]))

    # Parse the 6th line of the mineral species data
    words = f.readline().strip().split()
    data.Tmax = float(words[0])

    # Set the type of the species
    data.type = 'Mineral'

    return data

def correctAqueousData(data):
    # Correct how the charge of the species is represented in the formula
    data.formula = data.formula \
        .replace('(+0)', '').replace('(-0)', '').replace('(0)', '') \
        .replace('(+)', '+').replace('(-)', '-') \
        .replace('(+1)', '+1').replace('(-1)', '-') .replace('(1)', '+') \
        .replace('(+2)', '+2').replace('(-2)', '-2').replace('(2)', '+2') \
        .replace('(+3)', '+3').replace('(-3)', '-3').replace('(3)', '+3') \
        .replace('(+4)', '+4').replace('(-4)', '-4').replace('(4)', '+4')

    # Correct the name convention of the aqueous species from reference PC2 in SUPCRT92
    if data.references == 'PC2':
        data.name = data.name \
            .replace('1-', '-').replace('2-', '-2') \
            .replace('3-', '-3').replace('4-', '-4') \
            .replace('1' , '+') # necessary to change d+H3GDP1 to d+H3GDP+

    # Coorect the name convention of the aqueous species from reference PC3 in SUPCRT92 (keep the order below!)
    if data.references == 'PC3':
        data.name = data.name \
            .replace('1-red', 'red-').replace('1-ox', 'ox-') \
            .replace('2-red', 'red-2').replace('2-ox', 'ox-2') \
            .replace('3-red', 'red-3').replace('3-ox', 'ox-3') \
            .replace('4-red', 'red-4').replace('4-ox', 'ox-4') \
            .replace('-red', 'red-').replace('-ox', 'ox-') \
            .replace('+red', 'red+').replace('+ox', 'ox+')

        data.name = data.name \
            .replace('1-', '-').replace('2-', '-2') \
            .replace('3-', '-3').replace('4-', '-4')

    # Remove the following trailing suffixes from the name of the aqueous species
    data.name = data.name \
        .replace(',AQ', '') \
        .replace(',Aq', '') \
        .replace(',aq', '')

    # Correct the charge suffix of the aqueous charged species (keep the order below!)
    data.name = data.name[:-2] + data.name[-2:]

    # Ensure neutral aqueous species has the suffix (aq) and charged species has its suffix charge
    if data.charge == 0:
        data.name = data.name + '(aq)'
    if data.charge == 1 and data.name[-1:] != '+':
        data.name = data.name + '+'
    if data.charge == 2 and data.name[-2:] != '+2':
        data.name = data.name + '+2'
    if data.charge == 3 and data.name[-2:] != '+3':
        data.name = data.name + '+3'
    if data.charge == 4 and data.name[-2:] != '+4':
        data.name = data.name + '+4'
    if data.charge == -1 and data.name[-1:] != '-':
        data.name = data.name + '-'
    if data.charge == -2 and data.name[-2:] != '-2':
        data.name = data.name + '-2'
    if data.charge == -3 and data.name[-2:] != '-3':
        data.name = data.name + '-3'
    if data.charge == -4 and data.name[-2:] != '-4':
        data.name = data.name + '-4'

    # Correct the prefix of some aqueous species
    if data.name[:2] == 'A-':
        data.name = 'a-' + data.name[2:]
    if data.name[:2] == 'N-':
        data.name = 'n-' + data.name[2:]
    if data.name[:2] == 'O-':
        data.name = 'o-' + data.name[2:]
    if data.name[:2] == 'M-':
        data.name = 'm-' + data.name[2:]
    if data.name[:2] == 'P-':
        data.name = 'p-' + data.name[2:]

def correctGaseousData(data):
    # Correct the names of the gases
    data.gas = data.formula.title()
    data.name = data.abbreviation
    data.formula = data.abbreviation.replace('(g)', '')

def correctMineralData(data):
    # Correct the suffixes of some minerals
    data.name = data.name.replace(',High', ',high').replace(',Low', ',low')
    data.name = data.name.replace('-Ord', ',ord').replace('-Dis', ',dis').replace(',Ordered', ',ord')
    data.name = data.name.replace(',Dehydrated', ',dehydrated')
    data.name = data.name.replace(',Native', ',native')
    data.name = data.name.replace(',Hydrous', ',hydrous')
    data.name = data.name.replace(',Alpha', ',alpha')
    data.name = data.name.replace(',Beta', ',beta')
    data.name = data.name.replace(',Maximum', ',maximum')

    # Remove the symbol (s) from the name and formula of some minerals
    data.name = data.name.replace('(s)', '')
    data.formula = data.formula.replace('(s)', '')

    # Correct the thermodynamic data of the mineral
    data.Gf = None if data.Gf == 999999.0 else data.Gf
    data.Hf = None if data.Hf == 999999.0 else data.Hf

    if data.nptrans > 0:
        for i in range(data.nptrans):
            data.Htr[i]    = None if data.Htr[i] == 999999.0 else data.Htr[i]
            data.Vtr[i]    = None if data.Vtr[i] == 999999.0 else data.Vtr[i]
            data.dPdTtr[i] = None if data.dPdTtr[i] == 999999.0 else data.dPdTtr[i]

# Conversion constants
calorieToJoule = 4.184
barToPascal    = 1.0e+05
cm3Tom3        = 1.0e-06

def convertUnitsAqueousData(data):
    data.Gf *= calorieToJoule             # from cal/mol       to J/mol
    data.Hf *= calorieToJoule             # from cal/mol       to J/mol
    data.Sr *= calorieToJoule             # from cal/mol/K     to J/mol/K
    data.a1 *= calorieToJoule/barToPascal # from cal/mol/bar   to J/mol/Pa
    data.a2 *= calorieToJoule             # from cal/mol       to J/mol
    data.a3 *= calorieToJoule/barToPascal # from cal*K/mol/bar to J*K/mol/Pa
    data.a4 *= calorieToJoule             # from cal*K/mol     to J*K/mol
    data.c1 *= calorieToJoule             # from cal/mol/K     to J/mol/K
    data.c2 *= calorieToJoule             # from cal*K/mol     to J*K/mol
    data.wref *= calorieToJoule           # from cal/mol       to J/mol

def convertUnitsGaseousData(data):
    data.Gf *= calorieToJoule             # from cal/mol     to J/mol
    data.Hf *= calorieToJoule             # from cal/mol     to J/mol
    data.Sr *= calorieToJoule             # from cal/mol/K   to J/mol/K
    data.a  *= calorieToJoule             # from cal/mol/K   to J/mol/K
    data.b  *= calorieToJoule             # from cal/mol/K2 to J/mol/K2
    data.c  *= calorieToJoule             # from cal*K/mol   to J*K/mol

def convertUnitsMineralData(data):
    if data.Gf != None: data.Gf *= calorieToJoule # from cal/mol   to J/mol
    if data.Hf != None: data.Hf *= calorieToJoule # from cal/mol   to J/mol
    data.Sr *= calorieToJoule             # from cal/mol/K to J/mol/K
    data.Vr *= cm3Tom3                    # from cm3/mol to m3/mol
    if data.nptrans == 0:
        data.a *= calorieToJoule          # from cal/mol/K   to J/mol/K
        data.b *= calorieToJoule          # from cal/mol/K2 to J/mol/K2
        data.c *= calorieToJoule          # from cal*K/mol   to J*K/mol
    else:
        for i in range(data.nptrans):
            data.a[i] *= calorieToJoule   # from cal/mol/K   to J/mol/K
            data.b[i] *= calorieToJoule   # from cal/mol/K2 to J/mol/K2
            data.c[i] *= calorieToJoule   # from cal*K/mol   to J*K/mol

            if i < data.nptrans:
                if data.Htr[i] != None: data.Htr[i] *= calorieToJoule    # from cal/mol  to J/mol
                if data.Vtr[i] != None: data.Vtr[i] *= cm3Tom3           # from cm3/mol to m3/mol
                if data.dPdTtr[i] != None: data.dPdTtr[i] *= barToPascal # from bar/K    to Pa/K

def parseDatabase(filename):
    # Open the SUPCRT92 conforming database
    f = open(filename)

    # Create a list to store all the species data
    database = []

    # Parse the mineral species that do not undergo phase transitions
    jumpToSection(f, "minerals that do not undergo phase transitions")
    while True:
        data = parseMineralData(f, 0)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database.append(data)

    # Parse the mspecies species that undergo one phase transition
    jumpToSection(f, "minerals that undergo one phase transition");
    while True:
        data = parseMineralData(f, 1)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database.append(data)

    # Parse the mspecies species that undergo two phase transitions
    jumpToSection(f, "minerals that undergo two phase transitions");
    while True:
        data = parseMineralData(f, 2)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database.append(data)

    # Parse the mspecies species that undergo three phase transitions
    jumpToSection(f, "minerals that undergo three phase transitions");
    while True:
        data = parseMineralData(f, 3)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database.append(data)

    # Parse the gaseous species
    jumpToSection(f, "gases");
    while True:
        data = parseGaseousData(f)
        if data == None: break
        correctGaseousData(data)
##        convertUnitsGaseousData(data)
        database.append(data)

    # Parse the aspecies species
    jumpToSection(f, "aqueous species");
    while True:
        data = parseAqueousData(f)
        if data == None: break
        correctAqueousData(data)
##        convertUnitsAqueousData(data)
        data.dissociation = aqueous_complexes.get(data.name, '')
        database.append(data)
    return database


def createValueUnitsDict(value, units):
    res = OrderedDict()
    res['value'] = value
    res['units'] = units


def createGeneralSpecies(data):
    elemental_formula = ''.join([x[0] + '(' + str(x[1]) + ')' for x in data.elemental_formula])
    molar_mass = molarMass(data.elemental_formula)
    species = OrderedDict()
    species['Name'] = data.name
    species['Formula'] = data.formula
    species['Elements'] = elemental_formula
    species['Type'] = data.type
    species['MolarMass'] = createValueUnitsDict(molar_mass, 'g/mol')
    return species


def createWaterSpecies():
    species = OrderedDict()
    species['Name'] = 'H2O(l)'
    species['Formula'] = 'H2O'
    species['Elements'] = 'H(2)O(1)'
    species['Type'] = 'Aqueous'
    species['MolarMass'] = OrderedDict()
    species['MolarMass']['value'] = 18.0153
    species['MolarMass']['units'] = 'g/mol'
    species['Charge'] = 0
    return species


def createAqueousSpecies(data):
    species = OrderedDict()
    species['Charge'] = data.charge
    if data.dissociation != '':
        species['Dissociation'] = data.dissociation

    hkf = OrderedDict()
    hkf['References'] = data.references
    hkf['Date'] = data.date
    hkf['Gf'] = createValueUnitsDict(data.Gf, 'cal/mol')
    hkf['Hf'] = createValueUnitsDict(data.Hf, 'cal/mol')
    hkf['Sr'] = createValueUnitsDict(data.Sr, 'cal/(mol*K)')
    hkf['a1'] = createValueUnitsDict(data.a1, 'cal/(mol*bar)')
    hkf['a2'] = createValueUnitsDict(data.a2, 'cal/mol')
    hkf['a3'] = createValueUnitsDict(data.a3, '(cal*K)/(mol*bar)')
    hkf['a4'] = createValueUnitsDict(data.a4, '(cal*K)/mol')
    hkf['c1'] = createValueUnitsDict(data.c1, 'cal/(mol*K)')
    hkf['c2'] = createValueUnitsDict(data.c2, '(cal*K)/mol')
    hkf['wref'] = createValueUnitsDict(data.wref, 'cal/mol')

    species['Thermo'] = OrderedDict()
    species['Thermo']['HKF'] = hkf
    return species


def createGaseousSpecies(data):
    species = OrderedDict()
    species['GasName'] = data.gas

    # Get the critical properties of the current gaseous species
    properties = complement.critical_properties.get(data.name)
    if properties is not None:
        species['CriticalTemperature'] = createValueUnitsDict(properties['Tc']), 'K')
        species['CriticalPressure'] = createValueUnitsDict(properties['Pc']), 'bar')
        species['AcentricFactor'] = properties['omega']

    hkf = OrderedDict()
    hkf['References'] = data.references
    hkf['Date'] = data.date
    hkf['Gf'] = createValueUnitsDict(data.Gf, 'cal/mol')
    hkf['Hf'] = createValueUnitsDict(data.Hf, 'cal/mol')
    hkf['Sr'] = createValueUnitsDict(data.Sr, 'cal/(mol*K)')
    hkf['a'] = createValueUnitsDict(data.a, 'cal/(mol*K)')
    hkf['b'] = createValueUnitsDict(data.b, 'cal/(mol*K2)')
    hkf['c'] = createValueUnitsDict(data.c, '(cal*K)/mol')
    hkf['Tmax'] = createValueUnitsDict(data.Tmax, 'K')

    species['Thermo'] = OrderedDict()
    species['Thermo']['HKF'] = hkf
    return species


def createMineralSpecies(data):
    species = OrderedDict()
    species['MolarVolume'] = createValueUnitsDict(data.Vr, 'm3/mol')

    hkf = doc.createElement('hkf')
    appendElement(doc, hkf, 'references', str(data.references))
    appendElement(doc, hkf, 'date', str(data.date))
    appendElement(doc, hkf, 'Gf', str(data.Gf) if data.Gf != None else '', ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Hf', str(data.Hf) if data.Hf != None else '', ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Sr', str(data.Sr) if data.Sr != None else '', ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'Vr', str(data.Vr), ('units', 'cm3/mol'))
    appendElement(doc, hkf, 'nptrans', str(data.nptrans))
    appendElement(doc, hkf, 'Tmax', str(data.Tmax), ('units', 'K'))

    if data.nptrans == 0:
        appendElement(doc, hkf, 'a', str(data.a), ('units', 'cal/(mol*K)'))
        appendElement(doc, hkf, 'b', str(data.b), ('units', 'cal/(mol*K2)'))
        appendElement(doc, hkf, 'c', str(data.c), ('units', '(cal*K)/mol'))
    else:
        for i in range(data.nptrans+1):
            # Create a temperature range element for the thermodynamic data
            temperature_range = doc.createElement('temperature_range' + str(i))
            hkf.appendChild(temperature_range)
            appendElement(doc, temperature_range, 'a', str(data.a[i]), ('units', 'cal/(mol*K)'))
            appendElement(doc, temperature_range, 'b', str(data.b[i]), ('units', 'cal/(mol*K2)'))
            appendElement(doc, temperature_range, 'c', str(data.c[i]), ('units', '(cal*K)/mol'))

            if i < data.nptrans:
                appendElement(doc, temperature_range, 'Ttr', str(data.Ttr[i]))
                appendElement(doc, temperature_range, 'Htr', str(data.Htr[i]) if data.Htr[i] != None else '')
                appendElement(doc, temperature_range, 'Vtr', str(data.Vtr[i]) if data.Vtr[i] != None else '')
                appendElement(doc, temperature_range, 'dPdTtr', str(data.dPdTtr[i]) if data.dPdTtr[i] != None else '')

    thermo = doc.createElement('thermo')
    thermo.appendChild(hkf)
    species.appendChild(thermo)

# Parse the SUPCRT92 database
datalist = parseDatabase('slop98.reaktoro.dat')

# Collect the aqueous, gaseous and mineral species data
aqueous_datalist = [data for data in datalist if data.type == 'Aqueous']
gaseous_datalist = [data for data in datalist if data.type == 'Gaseous']
mineral_datalist = [data for data in datalist if data.type == 'Mineral']

if len(aqueous_datalist) + len(gaseous_datalist) + len(mineral_datalist) != len(datalist):
    raise RuntimeError('the number of aqueous, gaseous and mineral species does not match with the total number of species')

# Create the OrderedDict instance representing the whole document
doc = OrderedDict()

doc['Elements'] = []
doc['Species'] = []

# Insert the chemical elements in the OrderedDict `doc`
for (name, value) in sorted(complement.elements.items(), key=lambda x: x[1]):
    element = OrderedDict()
    element['Name'] = name
    element['MolarMass'] = OrderedDict()
    element['MolarMass']['value'] = value
    element['MolarMass']['units'] = 'g/mol'
    doc['Elements'].append(element)

# Insert the mineral species in the OrderedDict `doc`
for data in mineral_datalist:
    species = createMineralSpecies(data)
    doc['Species'].append(species)

# Insert the gaseous species in the OrderedDict `doc`
for data in gaseous_datalist:
    species = createGaseousSpecies(data)
    doc['Species'].append(species)

# Insert the water species in the OrderedDict `doc`
water = createWaterSpecies()
doc['Species'].append(water)

# Insert the aqueous species in the OrderedDict `doc`
for data in aqueous_datalist:
    species = createAqueousSpecies(data)
    doc['Species'].append(species)


# Output the database in XML format
f = open('supcrt98.xml', 'w')
f.write(doc.toprettyxml(indent="  "))