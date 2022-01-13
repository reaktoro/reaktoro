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
//   •
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
    MineralPhases minerals("Lim Cal Portlandite Gbs Gp Corundum Periclase Fe2O3 Na2O K2O "
                           "monocarbonate hemicarbonate hydrotalcite");
    // Note: C4AsH14, C4Ac0.5H12, C4AcH11 are not found
    
    // Define solid phases
    SolidPhase solidphase_C3AFS084H("C3FS0.84H4.32 C3AFS0.84H4.32");
    SolidPhase solidphase_ettringite_Al("ettringite Fe-ettringite");
    SolidPhase solidphase_monosulphate_Fe("Fe-monosulph05 Fe-monosulphate");
    SolidPhase solidphase_ettringite("ettringite ettringite30");
    SolidPhase solidphase_OH_SO4_AFm("C4AH19 monosulphate14");
    SolidPhase solidphase_CO3_SO4_AFt("tricarboalu03 ettringite03_ss");
    SolidPhase solidphase_CSHQ("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH");

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals,
                          solidphase_C3AFS084H,
                          solidphase_ettringite_Al,
                          solidphase_monosulphate_Fe,
                          solidphase_ettringite,
                          solidphase_OH_SO4_AFm,
                          solidphase_CSHQ);
    // Note: adding phase `solidphase_CO3_SO4_AFt` causes error Assertion failed: (values.minCoeff() >= 0.0), function setSpeciesAmounts while calling `solver.solve(state, conditions);` below

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    // Define equilibrium solver
    EquilibriumSolver solver(specs);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setSpeciesMass("H2O@",     63.9, "g");
    state.setSpeciesMass("Corundum",  4.9, "g"); // Al2O3
    state.setSpeciesMass("Lim",      63.9, "g"); // CaO
    state.setSpeciesMass("Fe2O3",     3.2, "g");
    state.setSpeciesMass("K2O",      0.78, "g");
    state.setSpeciesMass("Periclase", 1.8, "g"); // MgO
    state.setSpeciesMass("Na2O",     0.42, "g");
    state.setSpeciesMass("O2@",     2.145, "g"); // 1.0 (provided for O2) + 0.5*2.29 (provided for SO3)
    state.setSpeciesMass("SO2@",     2.29, "g");
    state.setSpeciesMass("SiO2@",    20.1, "g");
    state.setSpeciesMass("Ca(CO3)@",  0.1, "g");

    // Define temperature and pressure
    double T = 20.0; // in Celsius
    double P = 1.0; // in bar

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");

    // Equilibrate the initial state with given conditions and component amounts
    solver.solve(state, conditions);

    // Output the chemical state to a console
    state.output("state-cemdata.txt");

    return 0;
}