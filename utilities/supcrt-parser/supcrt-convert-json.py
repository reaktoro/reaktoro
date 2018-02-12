#!/usr/bin/env python

import argparse, sys
from xml.dom.minidom import Document
import json
from complement import *
from collections import OrderedDict

# The list of names writen in capitals that should be converted to tilte style
capital_names = set(open('capital-names.txt').read().splitlines())

def molarMass(elemental_formula):
    return sum([atoms * elements[name] for name, atoms in elemental_formula])

def parseElementalFormula(elemental_formula):
    # Remove the trailing symbols (aq) and (s) from the elemental formula of some aqueous and mineral species
    elemental_formula = elemental_formula.replace('(aq)', '')
    elemental_formula = elemental_formula.replace('(s)', '')

    # Split the elemental formula in words delimited by ( or )
    words = elemental_formula.replace(')', '(').split('(')
    words = [x for x in words if x != '']

    # Check if the elemental formula contains only one element
    if len(words) == 1:
        return [(words[0], 1)]

    # Collect the pairs (element, natoms) in the list elements
    elements = OrderedDict()
    for i in range(0, len(words), 2):
        element = words[i]
        natoms = int(words[i + 1])
        if element in ['+', '-'] and natoms != 0:
            elements['Z'] = natoms if element == '+' else -natoms
        else:
            elements[element] = natoms
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

    data['original_name'] = words[0] # the original and unchanged name
    data['name'] = words[0]
    data['formula'] = words[1]

    # Change the name of the species to title style if it is writen in caps
    if words[0] in capital_names:
        data['name'] = words[0].title()

    # Parse the 2nd line of the species data
    words = f.readline().strip().split()

    data['abbreviation'] = words[0]
    data['elements'] = parseElementalFormula(words[1])

    # Parse the 3rd line of the species data
    words = f.readline().strip().split()
    data['references'] = words[0].replace('ref:', '').replace('REF:', '')
    data['date'] = words[1].title()
    return True

def parseAqueousData(f):
    # Create the species data (see Table 11 of SUPCRT92 paper)
    data = OrderedDict()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Set the type of the species
    data['type'] = 'Aqueous'

    # Create a thermoProperties key
    thermo = OrderedDict()

    # Create a HKF key inside thermoProperties
    hkf = thermo['hkf'] = OrderedDict()

    # Parse the 4th line of the aqueous species data
    words = f.readline().strip().split()
    hkf['Gf'] = float(words[0])
    hkf['Hf'] = float(words[1])
    hkf['Sr'] = float(words[2])

    # Parse the 5th line of the aqueous species data
    words = f.readline().strip().split()
    hkf['a1'] = float(str(float(words[0]) * 1.0e-01))
    hkf['a2'] = float(str(float(words[1]) * 1.0e+02))
    hkf['a3'] = float(str(float(words[2]) * 1.0e+00))
    hkf['a4'] = float(str(float(words[3]) * 1.0e+04))

    # Parse the 6th line of the aqueous species data
    words = f.readline().strip().split()
    hkf['c1']     = float(str(float(words[0]) * 1.0e+00))
    hkf['c2']     = float(str(float(words[1]) * 1.0e+04))
    hkf['wref']   = float(str(float(words[2]) * 1.0e+05))
    data['charge'] = int(float(words[3]))

    # Set the thermoProperties entry of the species
    data['thermoProperties'] = thermo

    return data

def parseGaseousData(f):
    # Create the species data (see Table 9 of SUPCRT92 paper)
    data = OrderedDict()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Set the type of the species
    data['type'] = 'Gaseous'

    # Create a thermoProperties key
    thermo = OrderedDict()

    # Create a Maier-Kelly key inside thermoProperties
    maierkelly = thermo['maierkelly'] = OrderedDict()

    # Parse the 4th line of the gaseous species data
    words = f.readline().strip().split()
    maierkelly['Gf'] = float(words[0])
    maierkelly['Hf'] = float(words[1])
    maierkelly['Sr'] = float(words[2])

    # Parse the 5th line of the gaseous species data
    words = f.readline().strip().split()
    maierkelly['a'] = float(str(float(words[0]) * 1.0e+00))
    maierkelly['b'] = float(str(float(words[1]) * 1.0e-03))
    maierkelly['c'] = float(str(float(words[2]) * 1.0e+05))

    # Parse the 6th line of the gaseous species data
    words = f.readline().strip().split()
    maierkelly['Tmax'] = float(words[0])

    # Set the thermoProperties entry of the species
    data['thermoProperties'] = thermo

    return data

def parseMineralData(f, numPhaseTransitions):
    # Create the species data (see Table 10 of SUPCRT92 paper)
    data = OrderedDict()

    # Parse the common species data and check for continuation
    if not parseGeneralData(f, data):
        return None

    # Set the type of the species
    data['type'] = 'Mineral'

    # Create a thermoProperties key
    thermo = OrderedDict()

    # Create a Maier-Kelly key inside thermoProperties
    maierkelly = thermo['maierkelly'] = OrderedDict()

    # Parse the 4th line of the mineral species data
    words = f.readline().strip().split()
    maierkelly['Gf'] = float(words[0])
    maierkelly['Hf'] = float(words[1])
    maierkelly['Sr'] = float(words[2])
    maierkelly['Vr'] = float(words[3])
    maierkelly['numPhaseTransitions'] = numPhaseTransitions

    # Parse the 5th line of the mineral species data
    if numPhaseTransitions == 0:
        words = f.readline().strip().split()
        maierkelly['a'] = float(str(float(words[0]) * 1.0e+00))
        maierkelly['b'] = float(str(float(words[1]) * 1.0e-03))
        maierkelly['c'] = float(str(float(words[2]) * 1.0e+05))
    else:
        maierkelly['a'], maierkelly['b'], maierkelly['c'] = [], [], []
        maierkelly['Ttr'], maierkelly['Htr'], maierkelly['Vtr'], maierkelly['dPdTtr'] = [], [], [], []
        for i in range(numPhaseTransitions+1):
            words = f.readline().strip().split()
            maierkelly['a'].append(float(str(float(words[0]) * 1.0e+00)))
            maierkelly['b'].append(float(str(float(words[1]) * 1.0e-03)))
            maierkelly['c'].append(float(str(float(words[2]) * 1.0e+05)))

            if i < numPhaseTransitions:
                maierkelly['Ttr'].append(float(words[3]))
                maierkelly['Htr'].append(float(words[4]))
                maierkelly['Vtr'].append(float(words[5]))
                maierkelly['dPdTtr'].append(float(words[6]))

    # Parse the 6th line of the mineral species data
    words = f.readline().strip().split()
    maierkelly['Tmax'] = float(words[0])

    # Set the thermoProperties entry of the species
    data['thermoProperties'] = thermo

    return data

def correctAqueousData(data):
    # Correct the +1 and -1 suffixes in some species
    if data['name'][-2:] == '+1': data['name'] = data['name'][:-2] + '+'
    if data['name'][-2:] == '-1': data['name'] = data['name'][:-2] + '-'

    # Correct how the charge of the species is represented in the formula
    data['formula'] = data['formula'] \
        .replace('(+0)', '') \
        .replace('(-0)', '') \
        .replace('(+)', '+') \
        .replace('(-)', '-') \
        .replace('(+1)', '+') \
        .replace('(+2)', '+2') \
        .replace('(+3)', '+3') \
        .replace('(+4)', '+4') \
        .replace('(-1)', '-') \
        .replace('(-2)', '-2') \
        .replace('(-3)', '-3') \
        .replace('(-4)', '-4') \
        .replace('(0)', '') \
        .replace('(1)', '+') \
        .replace('(2)', '+2') \
        .replace('(3)', '+3') \
        .replace('(4)', '+4') \
        .replace('+1', '+') \
        .replace('-1', '-')

    # # Correct the name convention of the aqueous species from reference PC2 in SUPCRT92
    # if data.references == 'PC2':
    #     data.name = data.name \
    #         .replace('1-', '-') \
    #         .replace('2-', '-2') \
    #         .replace('3-', '-3') \
    #         .replace('4-', '-4') \
    #         .replace('1' , '+') # necessary to change d+H3GDP1 to d+H3GDP+

    # Coorect the name convention of the aqueous species from reference PC3 in SUPCRT92 (keep the order below!)
    if data['references'] == 'PC3':
        data['name'] = data['name'] \
            .replace('1-red', 'red-') \
            .replace('2-red', 'red--') \
            .replace('3-red', 'red---') \
            .replace('4-red', 'red----') \
            .replace('-red', 'red-') \
            .replace('+red', 'red+') \
            .replace('red', 'red') \
            .replace('1-ox', 'ox-') \
            .replace('2-ox', 'ox--') \
            .replace('3-ox', 'ox---') \
            .replace('4-ox', 'ox----') \
            .replace('-ox', 'ox-') \
            .replace('+ox', 'ox+')

        data['name'] = data['name'] \
            .replace('1-', '-') \
            .replace('2-', '--') \
            .replace('3-', '---') \
            .replace('4-', '----')

    # Remove the following trailing suffixes from the name of the aqueous species
    data['name'] = data['name'] \
        .replace(',AQ', '') \
        .replace(',Aq', '') \
        .replace(',aq', '')

    # Correct the charge suffix of the aqueous charged species (keep the order below!)
    data['name'] = data['name'][:-2] + data['name'][-2:]

    # Ensure neutral aqueous species has the suffix (aq) and charged species has its suffix charge
    if data['charge'] == 0:
        data['name'] = data['name'] + '(aq)'

    if data['references'] != 'PC2': # Skip the correction for ions in reference PC2 (the original names are messed up with wrong charges)
        data['name'] = data['name'] \
            .replace('+1', '+') \
            .replace('+2', '++') \
            .replace('+3', '+++') \
            .replace('+4', '++++') \
            .replace('-1', '-') \
            .replace('-2', '--') \
            .replace('-3', '---') \
            .replace('-4', '----')
        if data['charge'] == 1 and data['name'][-1:] != '+':
            data['name'] = data['name'] + '+'
        if data['charge'] == 2 and data['name'][-2:] != '++':
            data['name'] = data['name'] + '++'
        if data['charge'] == 3 and data['name'][-3:] != '+++':
            data['name'] = data['name'] + '+++'
        if data['charge'] == 4 and data['name'][-4:] != '++++':
            data['name'] = data['name'] + '++++'
        if data['charge'] == -1 and data['name'][-1:] != '-':
            data['name'] = data['name'] + '-'
        if data['charge'] == -2 and data['name'][-2:] != '--':
            data['name'] = data['name'] + '--'
        if data['charge'] == -3 and data['name'][-3:] != '---':
            data['name'] = data['name'] + '---'
        if data['charge'] == -4 and data['name'][-4:] != '----':
            data['name'] = data['name'] + '----'

    # Correct the prefix of some aqueous species
    if data['name'][:2] == 'A-':
        data['name'] = 'a-' + data['name'][2:]
    if data['name'][:2] == 'N-':
        data['name'] = 'n-' + data['name'][2:]
    if data['name'][:2] == 'O-':
        data['name'] = 'o-' + data['name'][2:]
    if data['name'][:2] == 'M-':
        data['name'] = 'm-' + data['name'][2:]
    if data['name'][:2] == 'P-':
        data['name'] = 'p-' + data['name'][2:]

def correctGaseousData(data):
    # Correct the names of the gases
    data['gas'] = data['formula'].title()
    data['name'] = data['abbreviation']
    data['formula'] = data['abbreviation'].replace('(g)', '')

    # Fix the names of some gaseous species
    if data['gas'] == 'Ortho-Cresol':
        data['name'] = 'o-Cresol(g)'
    if data['gas'] == 'Meta-Cresol':
        data['name'] = 'm-Cresol(g)'
    if data['gas'] == 'Para-Cresol':
        data['name'] = 'p-Cresol(g)'
    if data['gas'] == 'C2H4':
        data['gas'] = 'Ethylene'

def correctMineralData(data):
    # Correct the suffixes of some minerals
    data['name'] = data['name'].replace(',High', ',high').replace(',Low', ',low')
    data['name'] = data['name'].replace('-Ord', ',ord').replace('-Dis', ',dis').replace(',Ordered', ',ord').replace(',Disordered', ',dis')
    data['name'] = data['name'].replace(',Dehydrated', ',dehydrated')
    data['name'] = data['name'].replace(',Native', ',native')
    data['name'] = data['name'].replace(',Hydrous', ',hydrous')
    data['name'] = data['name'].replace(',Alpha', ',alpha')
    data['name'] = data['name'].replace(',Beta', ',beta')
    data['name'] = data['name'].replace(',Maximum', ',maximum')

    # Remove the symbol (s) from the name and formula of some minerals
    data['name'] = data['name'].replace('(s)', '')
    data['formula'] = data['formula'].replace('(s)', '')

    # Correct the thermodynamic data of the mineral
    maierkelly = data['thermoProperties']['maierkelly']

    if maierkelly['Gf'] == 999999.0: maierkelly['Gf'] = None
    if maierkelly['Hf'] == 999999.0: maierkelly['Hf'] = None

    if maierkelly['numPhaseTransitions'] > 0:
        for i in range(maierkelly['numPhaseTransitions']):
            if maierkelly['Htr'][i] == 999999.0: maierkelly['Htr'][i] = None
            if maierkelly['Vtr'][i] == 999999.0: maierkelly['Vtr'][i] = None
            if maierkelly['dPdTtr'][i] == 999999.0: maierkelly['dPdTtr'][i] = None

# Conversion constants
calorieToJoule = 4.184
barToPascal    = 1.0e+05
cm3Tom3        = 1.0e-06

def convertUnitsAqueousData(data):
    data['Gf'] *= calorieToJoule             # from cal/mol       to J/mol
    data['Hf'] *= calorieToJoule             # from cal/mol       to J/mol
    data['Sr'] *= calorieToJoule             # from cal/mol/K     to J/mol/K
    data['a1'] *= calorieToJoule/barToPascal # from cal/mol/bar   to J/mol/Pa
    data['a2'] *= calorieToJoule             # from cal/mol       to J/mol
    data['a3'] *= calorieToJoule/barToPascal # from cal*K/mol/bar to J*K/mol/Pa
    data['a4'] *= calorieToJoule             # from cal*K/mol     to J*K/mol
    data['c1'] *= calorieToJoule             # from cal/mol/K     to J/mol/K
    data['c2'] *= calorieToJoule             # from cal*K/mol     to J*K/mol
    data['wref'] *= calorieToJoule           # from cal/mol       to J/mol

def convertUnitsGaseousData(data):
    data['Gf'] *= calorieToJoule             # from cal/mol     to J/mol
    data['Hf'] *= calorieToJoule             # from cal/mol     to J/mol
    data['Sr'] *= calorieToJoule             # from cal/mol/K   to J/mol/K
    data['a']  *= calorieToJoule             # from cal/mol/K   to J/mol/K
    data['b']  *= calorieToJoule             # from cal/mol/K^2 to J/mol/K2
    data['c']  *= calorieToJoule             # from cal*K/mol   to J*K/mol

def convertUnitsMineralData(data):
    if data['Gf'] != None: data['Gf'] *= calorieToJoule # from cal/mol   to J/mol
    if data['Hf'] != None: data['Hf'] *= calorieToJoule # from cal/mol   to J/mol
    data['Sr'] *= calorieToJoule             # from cal/mol/K to J/mol/K
    data['Vr'] *= cm3Tom3                    # from cm3/mol to m3/mol
    if data['numPhaseTransitions'] == 0:
        data['a'] *= calorieToJoule          # from cal/mol/K   to J/mol/K
        data['b'] *= calorieToJoule          # from cal/mol/K^2 to J/mol/K2
        data['c'] *= calorieToJoule          # from cal*K/mol   to J*K/mol
    else:
        for i in range(data['numPhaseTransitions']):
            data['a'][i] *= calorieToJoule   # from cal/mol/K   to J/mol/K
            data['b'][i] *= calorieToJoule   # from cal/mol/K^2 to J/mol/K2
            data['c'][i] *= calorieToJoule   # from cal*K/mol   to J*K/mol

            if i < data['numPhaseTransitions']:
                if data['Htr'][i] != None: data['Htr'][i] *= calorieToJoule    # from cal/mol  to J/mol
                if data['Vtr'][i] != None: data['Vtr'][i] *= cm3Tom3           # from cm3/mol to m3/mol
                if data['dPdTtr'][i] != None: data['dPdTtr'][i] *= barToPascal # from bar/K    to Pa/K

def parseDatabase(filename):
    # Open the SUPCRT92 conforming database
    f = open(filename)

    # Create a list to store all the species data
    database = OrderedDict()

    database['elements'] = []
    database['species'] = []

    # Parse the mineral species that do not undergo phase transitions
    jumpToSection(f, "minerals that do not undergo phase transitions")
    while True:
        data = parseMineralData(f, 0)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database['species'].append(data)

    # Parse the mspecies species that undergo one phase transition
    jumpToSection(f, "minerals that undergo one phase transition");
    while True:
        data = parseMineralData(f, 1)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database['species'].append(data)

    # Parse the mspecies species that undergo two phase transitions
    jumpToSection(f, "minerals that undergo two phase transitions");
    while True:
        data = parseMineralData(f, 2)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database['species'].append(data)

    # Parse the mspecies species that undergo three phase transitions
    jumpToSection(f, "minerals that undergo three phase transitions");
    while True:
        data = parseMineralData(f, 3)
        if data == None: break
        correctMineralData(data)
##        convertUnitsMineralData(data)
        database['species'].append(data)

    # Parse the gaseous species
    jumpToSection(f, "gases");
    while True:
        data = parseGaseousData(f)
        if data == None: break
        correctGaseousData(data)
##        convertUnitsGaseousData(data)
        database['species'].append(data)

    # Parse the aspecies species
    jumpToSection(f, "aqueous species");
    while True:
        data = parseAqueousData(f)
        if data == None: break
        correctAqueousData(data)
##        convertUnitsAqueousData(data)
        data['dissociation'] = aqueous_complexes.get(data['name'], '')
        database['species'].append(data)

    database['criticalProperties'] = OrderedDict(critical_properties)

    # Create a units entry that lists the units of thermodynamic parameters
    database['units'] = OrderedDict()
    
    # Create an entry for the units of the HKF thermodynamic parameters
    database['units']['hkf'] = OrderedDict()
    database['units']['hkf']['Gf'] = 'cal/mol'
    database['units']['hkf']['Hf'] = 'cal/mol'
    database['units']['hkf']['Sr'] = 'cal/(mol*K)'
    database['units']['hkf']['a1'] = 'cal/(mol*bar)'
    database['units']['hkf']['a2'] = 'cal/mol'
    database['units']['hkf']['a3'] = '(cal*K)/(mol*bar)'
    database['units']['hkf']['a4'] = '(cal*K)/mol'
    database['units']['hkf']['c1'] = 'cal/(mol*K)'
    database['units']['hkf']['c2'] = '(cal*K)/mol'
    database['units']['hkf']['wref'] = 'cal/mol'

    # Create an entry for the units of the Maier-Kelly thermodynamic parameters
    database['units']['maierkelly'] = OrderedDict()
    database['units']['maierkelly']['Gf'] = 'cal/mol'
    database['units']['maierkelly']['Hf'] = 'cal/mol'
    database['units']['maierkelly']['Sr'] = 'cal/(mol*K)'
    database['units']['maierkelly']['Vr'] = 'cm3/mol'
    database['units']['maierkelly']['a'] = 'cal/(mol*K)'
    database['units']['maierkelly']['b'] = 'cal/(mol*K^2)'
    database['units']['maierkelly']['c'] = '(cal*K)/mol'
    database['units']['maierkelly']['Tmax'] = 'K'

    return database

def createElement(doc, name, text):
    elem = doc.createElement(name)
    elem.appendChild(doc.createTextNode(text))
    return elem

def appendElement(doc, root, name, text, attribute=None):
    elem = doc.createElement(name)
    if attribute != None: elem.setAttribute(attribute[0], attribute[1])
    elem.appendChild(doc.createTextNode(text))
    root.appendChild(elem)

def createSpeciesXML(doc, root, data):
    elemental_formula = ''.join([x[0] + '(' + str(x[1]) + ')' for x in data.elemental_formula])
    species = doc.createElement('Species')
    appendElement(doc, species, 'Name', data.name)
    appendElement(doc, species, 'Formula', data.formula)
    appendElement(doc, species, 'Elements', elemental_formula)
    appendElement(doc, species, 'Type', data.type)
    appendElement(doc, species, 'MolarMass', str(molarMass(data.elemental_formula)), ('units', 'g/mol'))
    root.appendChild(species)

    return species

def writeWaterSpeciesXML(doc, root):
    species = doc.createElement('Species')
    appendElement(doc, species, 'Name', 'H2O(l)')
    appendElement(doc, species, 'Formula', 'H2O')
    appendElement(doc, species, 'Elements', 'H(2)O(1)')
    appendElement(doc, species, 'Type', 'Aqueous')
    appendElement(doc, species, 'MolarMass', str(18.0153), ('units', 'g/mol'))
    appendElement(doc, species, 'Charge', str(0))
    root.appendChild(species)

def writeAqueousSpeciesXML(doc, root, data):
    species = createSpeciesXML(doc, root, data)
    appendElement(doc, species, 'Charge', str(data.charge))
    if data.dissociation != '':
        appendElement(doc, species, 'Dissociation', data.dissociation)

    hkf = doc.createElement('HKF')
    appendElement(doc, hkf, 'References', str(data.references))
    appendElement(doc, hkf, 'Date', str(data.date))
    appendElement(doc, hkf, 'Gf', str(data.Gf), ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Hf', str(data.Hf), ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Sr', str(data.Sr), ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'a1', str(data.a1), ('units', 'cal/(mol*bar)'))
    appendElement(doc, hkf, 'a2', str(data.a2), ('units', 'cal/mol'))
    appendElement(doc, hkf, 'a3', str(data.a3), ('units', '(cal*K)/(mol*bar)'))
    appendElement(doc, hkf, 'a4', str(data.a4), ('units', '(cal*K)/mol'))
    appendElement(doc, hkf, 'c1', str(data.c1), ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'c2', str(data.c2), ('units', '(cal*K)/mol'))
    appendElement(doc, hkf, 'wref', str(data.wref), ('units', 'cal/mol'))

    thermo = doc.createElement('Thermo')
    thermo.appendChild(hkf)
    species.appendChild(thermo)

def writeGaseousSpeciesXML(doc, root, data):
    species = createSpeciesXML(doc, root, data)
    appendElement(doc, species, 'Gas', str(data.gas))

    # Get the critical properties of the current gaseous species
    props = critical_properties.get(data.name)

    # Write them to XML nodes
    if props is not None:
        appendElement(doc, species, 'CriticalTemperature', str(props['Tc']), ('units', 'K'))
        appendElement(doc, species, 'CriticalPressure', str(props['Pc']), ('units', 'bar'))
        appendElement(doc, species, 'AcentricFactor', str(props['omega']))

    hkf = doc.createElement('HKF')
    appendElement(doc, hkf, 'References', str(data.references))
    appendElement(doc, hkf, 'Date', str(data.date))
    appendElement(doc, hkf, 'Gf', str(data.Gf), ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Hf', str(data.Hf), ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Sr', str(data.Sr), ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'a', str(data.a), ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'b', str(data.b), ('units', 'cal/(mol*K^2)'))
    appendElement(doc, hkf, 'c', str(data.c), ('units', '(cal*K)/mol'))
    appendElement(doc, hkf, 'Tmax', str(data.Tmax), ('units', 'K'))

    thermo = doc.createElement('Thermo')
    thermo.appendChild(hkf)
    species.appendChild(thermo)

def writeMineralSpeciesXML(doc, root, data):
    species = createSpeciesXML(doc, root, data)
    appendElement(doc, species, 'MolarVolume', str(data.Vr), ('units', 'm3/mol'))

    hkf = doc.createElement('HKF')
    appendElement(doc, hkf, 'References', str(data.references))
    appendElement(doc, hkf, 'Date', str(data.date))
    appendElement(doc, hkf, 'Gf', str(data.Gf) if data.Gf != None else '', ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Hf', str(data.Hf) if data.Hf != None else '', ('units', 'cal/mol'))
    appendElement(doc, hkf, 'Sr', str(data.Sr) if data.Sr != None else '', ('units', 'cal/(mol*K)'))
    appendElement(doc, hkf, 'Vr', str(data.Vr), ('units', 'cm3/mol'))
    appendElement(doc, hkf, 'NumPhaseTrans', str(data.numPhaseTransitions))
    appendElement(doc, hkf, 'Tmax', str(data.Tmax), ('units', 'K'))

    if data.numPhaseTransitions == 0:
        appendElement(doc, hkf, 'a', str(data.a), ('units', 'cal/(mol*K)'))
        appendElement(doc, hkf, 'b', str(data.b), ('units', 'cal/(mol*K^2)'))
        appendElement(doc, hkf, 'c', str(data.c), ('units', '(cal*K)/mol'))
    else:
        for i in range(data.numPhaseTransitions+1):
            # Create a temperature range element for the thermodynamic data
            temperature_range = doc.createElement('TemperatureRange' + str(i))
            hkf.appendChild(temperature_range)
            appendElement(doc, temperature_range, 'a', str(data.a[i]), ('units', 'cal/(mol*K)'))
            appendElement(doc, temperature_range, 'b', str(data.b[i]), ('units', 'cal/(mol*K^2)'))
            appendElement(doc, temperature_range, 'c', str(data.c[i]), ('units', '(cal*K)/mol'))

            if i < data.numPhaseTransitions:
                appendElement(doc, temperature_range, 'Ttr', str(data.Ttr[i]))
                appendElement(doc, temperature_range, 'Htr', str(data.Htr[i]) if data.Htr[i] != None else '')
                appendElement(doc, temperature_range, 'Vtr', str(data.Vtr[i]) if data.Vtr[i] != None else '')
                appendElement(doc, temperature_range, 'dPdTtr', str(data.dPdTtr[i]) if data.dPdTtr[i] != None else '')

    thermo = doc.createElement('Thermo')
    thermo.appendChild(hkf)
    species.appendChild(thermo)


def main():
    # Create a command-line argument parser
    parser = argparse.ArgumentParser(prog='ReaktoroParserSUPCRT')

    # Add the input argument
    parser.add_argument('input', type=str, \
        help='the relative path to the SUPCRT database file')

    # Add the output argument
    parser.add_argument('output', type=str, \
        help='the relative path of the output file')

    # Add the debug option (optional)
    parser.add_argument('-e', '--exclude', type=str, nargs='?', \
        help='the relative path to a file containing the names of species \
            to be excluded from the final database')

    # Parse the command-line arguments (remove the first argument, which is the name of this file
    args = parser.parse_args(sys.argv[1:])

    # Extract the command-line arguments
    infile = args.input
    outfile = args.output

    # Create a list of species that should be excluded
    excluded = []
    if args.exclude is not None:
        excluded = [line.strip() for line in open(args.exclude, 'r')]

    # Parse the SUPCRT92 database
    datalist = parseDatabase(infile)

    for s in datalist['species']:
        elements = dict(s['elements'])
        if elements.get('Z', 0) != s.get('charge', 0):
            print '{} has charge {} but Z {}'.format(s['name'], s.get('charge', 0), elements.get('Z', 0))

    # # Collect the aqueous, gaseous and mineral species data
    # aqueous_datalist = [data for data in datalist if data.type == 'Aqueous']
    # gaseous_datalist = [data for data in datalist if data.type == 'Gaseous']
    # mineral_datalist = [data for data in datalist if data.type == 'Mineral']

    # if len(aqueous_datalist) + len(gaseous_datalist) + len(mineral_datalist) != len(datalist):
    #     raise RuntimeError('the number of aqueous, gaseous and mineral species does not match with the total number of species')

    # # Remove from the excluded species from the lists above
    # aqueous_datalist = [x for x in aqueous_datalist if x.original_name not in excluded]
    # gaseous_datalist = [x for x in gaseous_datalist if x.original_name not in excluded]
    # mineral_datalist = [x for x in mineral_datalist if x.original_name not in excluded]

    # Output the database in XML format
    f = open(outfile, 'w')
    f.write(json.dumps(datalist, indent=4, separators=(',', ': ')))

    print 'File', outfile, 'has been successfully written!'
    

    # # Create the minidom document
    # doc = Document()

    # # Create the <database> element
    # db = doc.createElement("Database")
    # doc.appendChild(db)

    # # Create the xml elements for the chemical element
    # for (name, value) in sorted(elements.items(), key=lambda x: x[1]):
    #     element_node = doc.createElement("Element")
    #     appendElement(doc, element_node, 'Name', name)
    #     appendElement(doc, element_node, 'MolarMass', str(value), ('units', 'g/mol'))
    #     db.appendChild(element_node)

    # # Create the xml elements for the mineral species
    # for data in mineral_datalist:
    #     writeMineralSpeciesXML(doc, db, data)

    # # Create the xml elements for the gaseous species
    # for data in gaseous_datalist:
    #     writeGaseousSpeciesXML(doc, db, data)

    # # Create the xml elements for the aqueous species
    # writeWaterSpeciesXML(doc, db)
    # for data in aqueous_datalist:
    #     writeAqueousSpeciesXML(doc, db, data)

    # # Output the database in XML format
    # f = open(outfile, 'w')
    # f.write(doc.toprettyxml(indent="  "))

    # print 'File', outfile, 'has been successfully written!'

if __name__ == '__main__':
    main()