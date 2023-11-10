# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
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
import pytest


# TODO Implement tests for the python bindings of component AqueousProps in AqueousProps[test].py
def testAqueousProps():
    pass


@pytest.mark.parametrize(
    "db, composition", (
        (SupcrtDatabase("supcrt07"), {
            "H2O(aq)": (1, "kg"),
            "H+": (6.3e-9, "mol"),
            "OH-": (10e-6, "mol"),
            "Na+": (0.6, "mol"),
            "Cl-": (0.59767, "mol"),
            "CO2(aq)": (8.3e-6, "mol"),
            "HCO3-": (1660e-6, "mol"),
            "CO3-2": (331e-6, "mol"),
        }),
        (SupcrtDatabase("supcrtbl"), {
            "H2O(aq)": (1, "kg"),
            "H+": (6.3e-9, "mol"),
            "OH-": (10e-6, "mol"),
            "Na+": (0.6, "mol"),
            "Cl-": (0.59767, "mol"),
            "CO2(aq)": (8.3e-6, "mol"),
            "HCO3-": (1660e-6, "mol"),
            "CO3-2": (331e-6, "mol"),
        }),
        (PhreeqcDatabase("phreeqc.dat"), {
            "H2O": (1, "kg"),
            "H+": (6.3e-9, "mol"),
            "OH-": (10e-6, "mol"),
            "Na+": (0.6, "mol"),
            "Cl-": (0.59767, "mol"),
            "CO2": (8.3e-6, "mol"),
            "HCO3-": (1660e-6, "mol"),
            "CO3-2": (331e-6, "mol"),
        }),
        (ThermoFunDatabase("aq17"), {
            "H2O@": (1, "kg"),
            "H+": (6.3e-9, "mol"),
            "OH-": (10e-6, "mol"),
            "Na+": (0.6, "mol"),
            "Cl-": (0.59767, "mol"),
            "CO2@": (8.3e-6, "mol"),
            "HCO3-": (1660e-6, "mol"),
            "CO3-2": (331e-6, "mol"),
        })
    )
)
def testWolfGladrowSystemProps(db, composition):
    """
    This test checks if the pH and Alkalinity calculated with Reaktoro matches with the literature.

    The paper [1]_ provides a simple system with expected results for both pH and Alkalinity. In
    the paper, a titration with HCl is performed for the system. Here, we are omitting the
    titration and using the results for HCl = 0.0 mmol (see fig 1 of [1]_).

    .. [1] Wolf-Gladrow et al. (2007). See: https://doi.org/10.1016/j.marchem.2007.01.006
    """
    aq_species = " ".join(list(composition.keys()))
    aq_solution = AqueousPhase(aq_species)
    system = ChemicalSystem(db, aq_solution)

    state = ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "atm")
    for species_name, amount_and_units in composition.items():
        species_amount, species_unit = amount_and_units
        state.set(species_name, species_amount, species_unit)

    aq_props = AqueousProps(state)
    # Within 0.1% of relative deviation compared with the paper
    assert aq_props.pH().val() == pytest.approx(8.2, rel=1e-3)
    # Within 2% of relative deviation compared with the paper
    assert aq_props.alkalinity().val() == pytest.approx(2.33e-3, rel=2e-2)


def testAlkalinity():
    """
    This test uses the exact species that are considered in the Alkalinity calculation
    and checks if the expected result of the Total Alkalinity expression in [1]_ is
    obtained.
    
    .. [1] Wolf-Gladrow et al. (2007). See: https://doi.org/10.1016/j.marchem.2007.01.006
    """
    db = SupcrtDatabase("supcrt07")
    water = "H2O(aq) H+ OH-"
    ions = "Na+ Mg+2 Ca+2 K+ Sr+2 Cl- Br- NO3-"
    tpo4_cluster = "H3PO4(aq) H2PO4- HPO4-2 PO4-3"
    tnh3_cluster = "NH3(aq) NH4+"
    tso4_cluster = "SO4-2 HSO4-"
    thf_cluster = "F- HF(aq)"
    thno2_cluster = "NO2- HNO2(aq)"
    species_names = " ".join(
        [water, ions, tpo4_cluster, tnh3_cluster, tso4_cluster, thf_cluster, thno2_cluster]
    )
    aq_solution = AqueousPhase(species_names)
    system = ChemicalSystem(db, aq_solution)
    
    # Setting a fictional chemical state with known Alkalinity given Wolf-Gladrow formula
    state = ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "atm")
    state.set("H2O(aq)", 1, "kg")
    all_solute_species = species_names.split()
    all_solute_species.remove("H2O(aq)")
    for solute in all_solute_species:
        state.set(solute, 1, "mol")

    aq_props = AqueousProps(state)
    props = ChemicalProps(state)
    volume_in_liter = props.volume().val() * 1e3
    total_amounts_of_solutes_for_alkalinity = 1 + 2 + 2 + 1 + 2  # Na+ Mg+2 Ca+2 K+ Sr+2 
    total_amounts_of_solutes_for_alkalinity += -1 - 1 - 1 - 4  # Cl- Br- NO3- TPO4
    total_amounts_of_solutes_for_alkalinity += 2 - 4 - 2 - 2  # TNH3 TSO4 THF THNO2
    expected_alkalinity = total_amounts_of_solutes_for_alkalinity / volume_in_liter
    assert aq_props.alkalinity().val() == pytest.approx(expected_alkalinity)
