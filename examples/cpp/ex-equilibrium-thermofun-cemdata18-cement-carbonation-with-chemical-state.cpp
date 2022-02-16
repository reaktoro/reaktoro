// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Svetlana Kyas (27 September August 2021)
//
// and since revised by:
//   • G. Dan Miron (28 January 2022)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Define Thermofun database
    ThermoFunDatabase db("cemdata18");

    // Define aqueous phase
    AqueousPhase solution(speciate("H O K Na S Si Ca Mg Al C Cl"));

    // Set up a and b parameters for the ionic species (KOH, b = 0.123, a = 3.67)
    ActivityModelDebyeHuckelParams params;
    params.aiondefault = 3.67;
    params.biondefault = 0.123;
    params.bneutraldefault = 0.123;
    solution.setActivityModel(ActivityModelDebyeHuckel(params));

    // Define minerals
    MineralPhases minerals("Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag "
                           "C2S C3A C3S C4AF Lim Gp Brc K2SO4 K2O Na2SO4 Na2O");
    // Calcite, Hydrotalcite, Portlandite, Hemicarbonate, Monocarbonate, SiO2(am), Ferrihydrite, Gibbsite, Magnetite

    // Define solid phases
    SolidPhase solidphase_C3AFS084H("C3FS0.84H4.32 C3AFS0.84H4.32"); // AlFeSi-hydrogarnet_ss
    SolidPhase solidphase_ettringite("ettringite ettringite30"); // Ettrignite_ss
    SolidPhase solidphase_OH_SO4_AFm("C4AH13 monosulphate12"); // Monosulfate_ss
    SolidPhase solidphase_CSHQ("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH"); // CSH_ss

    // Define chemical system by providing database, aqueous phase, and minerals
    //ChemicalSystem system(db, solution, minerals);
    ChemicalSystem system(db, solution, minerals,
                          solidphase_C3AFS084H,
                          solidphase_ettringite,
                          solidphase_OH_SO4_AFm,
                          solidphase_CSHQ);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    // Define equilibrium solver
    EquilibriumOptions opts;

    // Define equilibrium solver and its result
    EquilibriumSolver solver(specs);
    EquilibriumResult res;

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setSpeciesMass("H2O@", 40, "g"); // water // water/binder = 0.4, 40g water per 100g of cement clinker
    // clinker phases
    state.setSpeciesMass("C2S" ,  9.70, "g"); // belite
    state.setSpeciesMass("C3A" ,  7.72, "g"); // aluminate
    state.setSpeciesMass("C3S" , 67.31, "g"); // alite
    state.setSpeciesMass("C4AF",  8.14, "g"); // ferrite, (CaO)4(Al2O3)(Fe|3|2O3)
    // additional
    state.setSpeciesMass("Gp"    , 3.13, "g"); // gypsum, CaSO4(H2O)2
    state.setSpeciesMass("Cal"   , 0.10, "g"); // calcite, CaCO3
    state.setSpeciesMass("Lim"   , 0.93, "g"); // lime, CaO
    state.setSpeciesMass("Brc"   , 1.31, "g"); // brucite, Mg(OH)2
    state.setSpeciesMass("K2SO4" , 1.34, "g"); // potasium-sulfate
    state.setSpeciesMass("K2O"   , 0.05, "g"); // potasium oxide
    state.setSpeciesMass("Na2SO4", 0.21, "g"); // sodium sulfate
    state.setSpeciesMass("Na2O"  , 0.05, "g"); // sodium oxide
    state.setSpeciesMass("O2@"   , 0.15, "g"); // oxygen to stabilize the system

    // Define temperature and pressure
    double T = 20.0; // in Celsius
    double P = 1.0; // in bar

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(state, conditions);
    std::cout << "res (cemdata18) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a text-file
    state.output("state-cemdata18.txt");

    // Output the chemical properties to a text-file
    ChemicalProps props(state);
    props.output("props.txt");

    // Output the aqueous properties to a text-file
    AqueousProps aprops(state);
    aprops.output("aprops.txt");

    return 0;
}