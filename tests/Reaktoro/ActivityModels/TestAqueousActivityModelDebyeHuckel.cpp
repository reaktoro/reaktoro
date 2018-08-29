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

TEST_CASE("Electrolyte Solution: NaCl")
{
	const Database db("supcrt98");
	ChemicalEditor editor(db);
	editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl-")
		.setChemicalModelDebyeHuckel();

	ChemicalSystem system(editor);

	ChemicalState state(system);
	state.setSpeciesAmount("H2O(l)", 1.0 / waterMolarMass);
	state.setSpeciesAmount("H+", 1.0e-7);
	state.setSpeciesAmount("OH-", 1.0e-7);
	state.setSpeciesAmount("Na+", 0.5);
	state.setSpeciesAmount("Cl-", 0.5);

	ChemicalProperties properties = state.properties();

	const Vector ln_g = properties.lnActivityCoefficients().val;

	const double ln10 = std::log(10.0);
	const double I = 0.5;
	const double A = 0.51137763214615395;
	const double B = 0.3287812566783076;
	//    const double nwo = 55.508472036052972;

	std::map<std::string, double> a = {
		{"H+", 9.0}, {"OH-", 3.5}, {"Na+", 4.08}, {"Cl-", 3.63}
	};

	std::map<std::string, double> b = {
		{"H+", 0.0}, {"OH-", 0.0}, {"Na+",  0.082}, {"Cl-", 0.017}
	};

	Vector ln_g_expected(5);
	ln_g_expected[1] = ln10 * (-A * std::sqrt(I) / (1 + B * a["H+"] * std::sqrt(I)) + b["H+"] * I);
	ln_g_expected[2] = ln10 * (-A * std::sqrt(I) / (1 + B * a["OH-"] * std::sqrt(I)) + b["OH-"] * I);
	ln_g_expected[3] = ln10 * (-A * std::sqrt(I) / (1 + B * a["Na+"] * std::sqrt(I)) + b["Na+"] * I);
	ln_g_expected[4] = ln10 * (-A * std::sqrt(I) / (1 + B * a["Cl-"] * std::sqrt(I)) + b["Cl-"] * I);
	
	CHECK(ln_g[5] == approx(ln_g_expected[1]));
	CHECK(ln_g[17] == approx(ln_g_expected[2]));
	CHECK(ln_g[13] == approx(ln_g_expected[3]));
	CHECK(ln_g[0] == approx(ln_g_expected[4]));
}

//TEST_CASE("Electrolyte Solution: KCl")
//{
//    DebyeHuckelParams db;
//
//    ChemicalEditor editor(db);
//    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ HCO3- CO2(aq) CO3--")
//        .setChemicalModelDebyeHuckel();
//
//    ChemicalSystem system(editor);
//
//    ChemicalState state(system);
//    state.setSpeciesAmount("H2O(l)", 1.0/waterMolarMass);
//    state.setSpeciesAmount("H+", 1.0e-4);
//    state.setSpeciesAmount("OH-", 1.0e-10);
//    state.setSpeciesAmount("Na+", 0.5);
//    state.setSpeciesAmount("Cl-", 0.5);
//    state.setSpeciesAmount("HCO3-", 0.1001);
//    state.setSpeciesAmount("CO2(aq)", 0.4);
//    state.setSpeciesAmount("CO3--", 1.0e-6);
//
//    ChemicalProperties properties = state.properties();
//
//    const double I = 0.55010200005;
//    const double A = 0.51137763214615395;
//    const double B = 0.3287812566783076;
//    const double nwo = 55.508472036052972;
//
//    const std::map<std::string, double> aion = {
//        {"H+", 9.0}, {"OH-", 3.5}, {"Na+", 4.08}, {"Cl-", 3.63},
//    };
//
//    const std::map<std::string, double> bion = {
//        {"H+", 0.0}, {"OH-", 0.0}, {"Na+",  0.082}, {"Cl-", 0.017}, {"HCO3-", 5.4}
//    };
//
//}
//
//TEST_CASE("Electrolyte Solution: CaCl2")
//{
//    DebyeHuckelParams db;
//
//    ChemicalEditor editor(db);
//    editor.addAqueousPhase("H2O(l) H+ OH- Ca++ Cl-")
//        .setChemicalModelDebyeHuckel();
//
//    ChemicalSystem system(editor);
//
//    ChemicalState state(system);
//    state.setSpeciesAmount("H2O(l)", 1.0/waterMolarMass);
//    state.setSpeciesAmount("H+", 1.0e-7);
//    state.setSpeciesAmount("OH-", 1.0e-7);
//    state.setSpeciesAmount("Ca++", 0.1);
//    state.setSpeciesAmount("Cl-", 0.2);
//
//    ChemicalProperties properties = state.properties();
//
//    const double I = 0.5*(0.1*4 + 0.2*1);
//    const double A = 0.51137763214615395;
//    const double B = 0.3287812566783076;
//    const double nwo = 55.508472036052972;
//
//    const std::map<std::string, double> aion = {
//        {"H+", 9.0}, {"OH-", 3.5}, {"Na+", 4.08}, {"Cl-", 3.63},
//    };
//
//    const std::map<std::string, double> bion = {
//        {"H+", 0.0}, {"OH-", 0.0}, {"Na+",  0.082}, {"Cl-", 0.017}, {"HCO3-", 5.4}
//    };
//
//}
//
//TEST_CASE("Electrolyte Solution: CO2")
//{
//    DebyeHuckelParams db;
//
//    ChemicalEditor editor(db);
//    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- HCO3- CO2(aq) CO3--")
//        .setChemicalModelDebyeHuckel();
//
//    ChemicalSystem system(editor);
//
//    ChemicalState state(system);
//    state.setSpeciesAmount("H2O(l)", 1.0/waterMolarMass);
//    state.setSpeciesAmount("H+", 1.0e-4);
//    state.setSpeciesAmount("OH-", 1.0e-10);
//    state.setSpeciesAmount("Na+", 0.5);
//    state.setSpeciesAmount("Cl-", 0.5);
//    state.setSpeciesAmount("HCO3-", 0.1001);
//    state.setSpeciesAmount("CO2(aq)", 0.4);
//    state.setSpeciesAmount("CO3--", 1.0e-6);
//
//    ChemicalProperties properties = state.properties();
//
//    const double I = 0.55010200005;
//    const double A = 0.51137763214615395;
//    const double B = 0.3287812566783076;
//    const double nwo = 55.508472036052972;
//
//    const std::map<std::string, double> aion = {
//        {"H+", 9.0}, {"OH-", 3.5}, {"Na+", 4.08}, {"Cl-", 3.63},
//    };
//
//    const std::map<std::string, double> bion = {
//        {"H+", 0.0}, {"OH-", 0.0}, {"Na+",  0.082}, {"Cl-", 0.017}, {"HCO3-", 5.4}
//    };
//
//}