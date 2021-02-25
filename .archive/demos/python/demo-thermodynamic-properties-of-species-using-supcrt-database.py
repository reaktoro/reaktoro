# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *
from numpy import *

# Create a Database object loaded with SUPCRT98 database file.
database = Database("supcrt98.xml")

# Create a Thermo object for thermodynamic property calculations
thermo = Thermo(database)

# Create shorter alias for the methods in Thermo that calculate
# standard partial molar Gibbs energy and enthalpy of species
evalG0 = thermo.standardPartialMolarGibbsEnergy
evalH0 = thermo.standardPartialMolarEnthalpy

# Create a list of species names, as found in the database
species = ['H2O(l)', 'HCO3-', 'CO2(aq)', 'CO2(g)', 'Calcite']

# Create a numeric array of temperature values (in units of K)
temperatures = array([25, 50, 75, 100, 200, 300]) + 273.15

# Create pressure variable (in units of Pa)
P = 100.e+5

# Create Python dictionaries containing the standard partial
# molar Gibbs energy and enthalpy for each species
G0 = {}
H0 = {}

# Calculate the standard chemical potentials of the species at
# 100 bar and at each temperature in the array of temperatures
for name in species:
    G0[name] = array([evalG0(T, P, name).val for T in temperatures])
    H0[name] = array([evalH0(T, P, name).val for T in temperatures])

# For each species, create a file containing three columns:
# 1st column: temperature (in units of K)
# 2nd column: standard partial molar gibbs energy (in units of J/mol)
# 3rd column: standard partial molar enthalpy (in units of J/mol)
for name in species:
    f = open('calculated-standard-species-properties-' + name + '.txt', 'w')
    print('T(K), G0(J/mol), H0(J/mol)', file=f)
    for Tval, G0val, H0val in zip(temperatures, G0[name], H0[name]):
        print('{0}, {1}, {2}'.format(Tval, G0val, H0val), file=f)
