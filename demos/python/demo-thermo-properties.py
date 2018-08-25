# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
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
import matplotlib.pyplot as plt

# Create a Database object loaded with SUPCRT database file.
database = Database("supcrt98.xml")

# Create a Thermo object for thermodynamic property calculations
thermo = Thermo(database)

# Define two reactions for which we'll calculate their log(K)
reaction1 = 'CO2(g) + H2O(l) = HCO3- + H+'
reaction2 = 'Calcite + H+ = Ca++ + HCO3-'

# Define pressure as 100 bar for the calculations (converted to Pa)
P = 100.0e5

# Create an array with temperature values from 25 to 300 C (converted to K)
x = linspace(25.0, 300.0, 50) + 273.15

# Calculate the log(K) of the reactions at those temperature points and 100 bar
y1 = [thermo.logEquilibriumConstant(T, P, reaction1).val for T in x]
y2 = [thermo.logEquilibriumConstant(T, P, reaction2).val for T in x]

x = x - 273.15

# Plot the log(K) of the reactions
plt.xlim(xmin=25, xmax=300)
plt.xlabel('Temperature [C]')
plt.ylabel('log(K)')
plt.plot(x, y1, label=reaction1)
plt.plot(x, y2, label=reaction2)
plt.show()
