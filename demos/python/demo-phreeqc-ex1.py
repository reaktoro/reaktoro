# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2015 Allan Leal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

# This string defines a PHREEQC script problem.
# This problem was taken from the official PHREEQC example named ex1.
ex1 = r'''(
TITLE Example 1.--Add uranium and speciate seawater.
SOLUTION 1  SEAWATER FROM NORDSTROM AND OTHERS (1979)
        units   ppm
        pH      8.22
        pe      8.451
        density 1.023
        temp    25.0
        redox   O(0)/O(-2)
        Ca              412.3
        Mg              1291.8
        Na              10768.0
        K               399.1
        Fe              0.002
        Mn              0.0002  pe
        Si              4.28
        Cl              19353.0
        Alkalinity      141.682 as HCO3
        S(6)            2712.0
        N(5)            0.29    gfw   62.0
        N(-3)           0.03    as    NH4
        U               3.3     ppb   N(5)/N(-3)
        O(0)            1.0     O2(g) -0.7
SOLUTION_MASTER_SPECIES
        U       U+4     0.0     238.0290     238.0290
        U(4)    U+4     0.0     238.0290
        U(5)    UO2+    0.0     238.0290
        U(6)    UO2+2   0.0     238.0290
SOLUTION_SPECIES
        #primary master species for U
        #is also secondary master species for U(4)
        U+4 = U+4
                log_k          0.0
        U+4 + 4 H2O = U(OH)4 + 4 H+
                log_k          -8.538
                delta_h        24.760 kcal
        U+4 + 5 H2O = U(OH)5- + 5 H+
                log_k          -13.147
                delta_h        27.580 kcal
        #secondary master species for U(5)
        U+4 + 2 H2O = UO2+ + 4 H+ + e-
                log_k          -6.432
                delta_h        31.130 kcal
        #secondary master species for U(6)
        U+4 + 2 H2O = UO2+2 + 4 H+ + 2 e-
                log_k          -9.217
                delta_h        34.430 kcal
        UO2+2 + H2O = UO2OH+ + H+
                log_k          -5.782
                delta_h        11.015 kcal
        2UO2+2 + 2H2O = (UO2)2(OH)2+2 + 2H+
                log_k          -5.626
                delta_h        -36.04 kcal
        3UO2+2 + 5H2O = (UO2)3(OH)5+ + 5H+
                log_k          -15.641
                delta_h        -44.27 kcal
        UO2+2 + CO3-2 = UO2CO3
                log_k          10.064
                delta_h        0.84 kcal
        UO2+2 + 2CO3-2 = UO2(CO3)2-2
                log_k          16.977
                delta_h        3.48 kcal
        UO2+2 + 3CO3-2 = UO2(CO3)3-4
                log_k          21.397
                delta_h        -8.78 kcal
PHASES
        Uraninite
        UO2 + 4 H+ = U+4 + 2 H2O
        log_k          -3.490
        delta_h        -18.630 kcal
END
)'''

# Initialize a Phreeqc instance with the official phreeqc.dat database file
phreeqc = Phreeqc('../../databases/phreeqc/phreeqc.dat')

# Execute a PHREEQC script defining a geochemical problem.
# Here this script is actually embedded into a string named `ex1`.
# However, `ex1` could also be a string containing the path to a script file.
# Method execute will automatically identify when the contents are embedded in
# the string and when the string is actually a path to a script file.
phreeqc.execute(ex1)

# Initialize a ChemicalSystem instance using the current state of the Phreeqc
# instance. This will allow the use of both PHREEQC thermodynamic data and
# PHREEQC activity models in the subsequent equilibrium calculations using
# Reaktoro's algorithms.
system = ChemicalSystem(phreeqc)

# Initialize an ChemicalState instance using the current state of the
# Phreeqc instance.
state = phreeqc.state(system)

# Output the equilibrium state calculated by PHREEQC to a file.
state.output('state-phreeqc.txt')

# Define an equilibrium problem in which the current state is mixed with 1
# mmol of HCl.
problem = EquilibriumProblem(system)
problem.setTemperature(state.temperature())
problem.setPressure(state.pressure())
problem.add(state)
problem.add('HCl', 1.0, 'mmol')

# Set Hessian of Gibbs energy to an approximation, since PHREEQC does not 
# compute molar derivatives of activities 
options = EquilibriumOptions()
options.hessian = GibbsHessian.Approximation

# Calculate the new equilibrium state of the system.
# This will use both PHREEQC thermodynamic data and PHREEQC activity models.
state = ChemicalState(system)
equilibrate(state, problem, options)

# Print the new equilibrium state and check with pH is more acidic now.
state.output('state-phreeqc-updated.txt')


