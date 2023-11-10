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
