# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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

# -----------------------------------------------------------------------------
# 👏 Acknowledgements 👏
# -----------------------------------------------------------------------------
# This example was originally authored by:
#   • Allan Leal (15 November 2022)
#
# and since revised by:
#   • Allan Leal (28 August 2023)
#     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
# -----------------------------------------------------------------------------


from reaktoro import *
from reaktplot import *


# The surface area of a cube per volume (in m2/m3)
Abar = 6.0

# The dissolution rate model for halite (NaCl)
def ratefn(props: ChemicalProps):
    aprops = AqueousProps(props)
    k0 = pow(10.0, -0.21)  # reaction rate constant at 25 °C from Palandri and Kharaka (2004)
    q = props.phaseProps("Halite").volume()  # volume in m3
    Omega = aprops.saturationRatio("Halite")
    return q * Abar * k0 * (1 - Omega)

db = PhreeqcDatabase("phreeqc.dat")

system = ChemicalSystem(db,
    AqueousPhase("H2O H+ OH- Na+ Cl-").set(ActivityModelPhreeqc(db)),
    MineralPhase("Halite"),
    GeneralReaction("Halite = Na+ + Cl-").setRateModel(ratefn)
)

state = ChemicalState(system)
state.set("H2O", 1.0, "kg")
state.scalePhaseVolume("Halite", 1.0, "cm3")  # start with 1 cm3 cube of halite crystal

solver = KineticsSolver(system)

table = Table()

dt = 60.0  # time step (in seconds)

# Initiate the time stepping for the kinetics modeling
for i in range(501):
    result = solver.solve(state, dt)

    assert result.succeeded(), f"Calculation did not succeed at time step #{i}."

    props = state.props()

    table.column("Time")   << i*dt / 60  # from seconds to minutes
    table.column("Halite") << props.phaseProps("Halite").volume() * 1e+6  # from m3 to cm3
    table.column("Na+")    << props.speciesAmount("Na+")  # in mol
    table.column("Cl-")    << props.speciesAmount("Cl-")  # in mol

table.save("ex-kinetics-halite-using-custom-rate-model.txt")

fig1 = Figure()
fig1.title("AQUEOUS SPECIES AMOUNTS OVER TIME")
fig1.xaxisTitle("Time [minute]")
fig1.yaxisTitle("Amount [mol]")
fig1.drawLine(table["Time"], table["Na+"], "Na<sup>+</sup>")
fig1.drawLine(table["Time"], table["Cl-"], "Cl<sup>-</sup>")
fig1.show()
fig1.save("ex-kinetics-halite-using-custom-rate-model-fig1.pdf")

fig2 = Figure()
fig2.title("CALCITE VOLUME OVER TIME")
fig2.xaxisTitle("Time [minute]")
fig2.yaxisTitle("Volume [m<sup>3</sup>]")
fig2.drawLine(table["Time"], table["Halite"], "Halite")
fig2.show()
fig2.save("ex-kinetics-halite-using-custom-rate-model-fig2.pdf")
