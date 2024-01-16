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
#   • Svetlana Kyas (14 July 2021)
# -----------------------------------------------------------------------------


from reaktoro import *


def printSpeciesInPhase(phase):
    phase_species = ""
    for species in phase.species():
        phase_species += species.name() + " "
    print("Species in " + phase.name() + ": " +  phase_species[0:-1])


# Definition of custom database with manula addition of the species (and their name, tags, ect)
db = Database()
db.addSpecies( Species("H2O(aq)") )
db.addSpecies( Species("H+") )
db.addSpecies( Species("OH-") )
db.addSpecies( Species("H2(aq)") )
db.addSpecies( Species("O2(aq)") )
db.addSpecies( Species("Na+") )
db.addSpecies( Species("Cl-") )
db.addSpecies( Species("NaCl(aq)") )
db.addSpecies( Species("HCl(aq)") )
db.addSpecies( Species("NaOH(aq)") )
db.addSpecies( Species("Ca++") )
db.addSpecies( Species("Mg++") )
db.addSpecies( Species("CO2(aq)") )
db.addSpecies( Species("HCO3-") )
db.addSpecies( Species("CO3--") )
db.addSpecies( Species("CaCl2(aq)") )
db.addSpecies( Species("MgCl2(aq)") )
db.addSpecies( Species("SiO2(aq)") )
db.addSpecies( Species("NaCl(s)").withName("Halite") )
db.addSpecies( Species("SiO2(s)").withName("Quartz") )
db.addSpecies( Species("CO2(g)") )
db.addSpecies( Species("O2(g)") )
db.addSpecies( Species("H2(g)") )
db.addSpecies( Species("H2O(g)") )
db.addSpecies( Species("CH4(g)") )
db.addSpecies( Species("CO(g)") )
db.addSpecies( Species("CaCO3(s)").withName("Calcite").withTags("carbonate") )
db.addSpecies( Species("MgCO3(s)").withName("Magnesite").withTags("carbonate") )
db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite").withTags("carbonate") )
db.addSpecies( Species("C(s)").withName("Graphite") )
db.addSpecies( Species("CaO(s)").withName("Lime") )
db.addSpecies( Species("N2(g)").withTags("inert") )
db.addSpecies( Species("C4H9OH").withTags("organic").withName("1-Butanol(aq)") )
db.addSpecies( Species("C4H8").withTags("organic").withName("1-Butene(aq)") )
db.addSpecies( Species("BaSO4(s)").withName("Barite").withTags("sulfate") )
db.addSpecies( Species("SrSO4(s)").withName("Celestite").withTags("sulfate") )
db.addSpecies( Species("PbSO4(s)").withName("Anglesite").withTags("sulfate") )
db.addSpecies( Species("CaSO4(s)").withName("Anhydrite").withTags("sulfate") )

print("# -----------------------------------------------------------------------------------------------------------------")
print("# Illustrate how select all the aqueous species with provided elements from database")
print("# -----------------------------------------------------------------------------------------------------------------")

# Define phases containing aqueous and gaseous species
phases = Phases(db)
phases.add( AqueousPhase(speciate("H O C Na Cl")) )
phases.add( GaseousPhase("CO2(g)") )

# Print species contained in these phases
phases = phases.convert()
printSpeciesInPhase(phases[0])
printSpeciesInPhase(phases[1])

print("# -----------------------------------------------------------------------------------------------------------------")
print("# Illustrate how to exclude aqueous species with a given tag")
print("# -----------------------------------------------------------------------------------------------------------------")

phases = Phases(db)
phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
phases.add( GaseousPhase("CO2(g)") )

phases = phases.convert()
printSpeciesInPhase(phases[0])
printSpeciesInPhase(phases[1])

print("# -----------------------------------------------------------------------------------------------------------------")
print("# Illustrate how to exclude aqueous species and minerals with given tags")
print("# -----------------------------------------------------------------------------------------------------------------")

phases = Phases(db)
phases.add( AqueousPhase(speciate("H O"), exclude("organic")))
phases.add( MineralPhases(speciate("C Ca O"), exclude("carbonate")) )

phases = phases.convert()
printSpeciesInPhase(phases[0])
printSpeciesInPhase(phases[1])
printSpeciesInPhase(phases[2])

print("# -----------------------------------------------------------------------------------------------------------------")
print("# Illustrate how to exclude aqueous and gaseous species as well as minerals with given tags")
print("# -----------------------------------------------------------------------------------------------------------------")

phases = Phases(db)
phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
phases.add( GaseousPhase(exclude("inert")) )
phases.add( MineralPhases(exclude("carbonate")) )

phases = phases.convert()
printSpeciesInPhase(phases[0])
printSpeciesInPhase(phases[1])
printSpeciesInPhase(phases[2])
printSpeciesInPhase(phases[3])
