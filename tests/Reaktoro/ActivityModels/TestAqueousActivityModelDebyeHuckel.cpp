// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2017 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <doctest/doctest.hpp>

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

auto createAqueousMixture() -> AqueousPhase
{
    ChemicalEditor editor;
    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ HCO3- CO2(aq) CO3--")
        .setChemicalModelDebyeHuckel();

    ChemicalSystem system;
}

TEST_CASE("TestAqueousActivityModelDebyeHuckel")
{
    ChemicalEditor editor;
    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ HCO3- CO2(aq) CO3--")
        .setChemicalModelDebyeHuckel();

    ChemicalSystem system(editor);

    ChemicalState state(system);
    state.setSpeciesAmount("H2O(l)", 1.0/waterMolarMass);
    state.setSpeciesAmount("H+", 1.0e-4);
    state.setSpeciesAmount("OH-", 1.0e-10);
    state.setSpeciesAmount("Na+", 0.5);
    state.setSpeciesAmount("Cl-", 0.5);
    state.setSpeciesAmount("HCO3-", 0.1001);
    state.setSpeciesAmount("CO2(aq)", 0.4);
    state.setSpeciesAmount("CO3--", 1.0e-6);

    ChemicalProperties properties = state.properties();

}
