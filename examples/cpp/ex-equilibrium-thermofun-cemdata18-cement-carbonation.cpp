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
    AqueousPhase aq_solution(speciate("H O K Na S Si Ca Mg Al C Cl"));

    // Define gas phase
    GaseousPhase gaseous(speciate("H O C"));

    // Here we define the activity model used to calculate the activity coefficients of the aqueous species
    // Set up a and b parameters for KOH background electrolyte b = 0.123, a = 3.67
    ActivityModelDebyeHuckelParams params;
    params.aiondefault = 3.67;
    params.biondefault = 0.123;
    params.bneutraldefault = 0.123;
    aq_solution.setActivityModel(ActivityModelDebyeHuckel(params));

    // Define pure mineral phases
    // Calcite, Hydrotalcite, Portlandite, Hemicarbonate, Monocarbonate, SiO2(am), Ferrihydrite, Gibbsite, Magnetite
    MineralPhases pure_minerals("Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag");
    

    // Define ideal solid-solution phases
    SolidPhase ss_C3AFS084H("C3FS0.84H4.32 C3AFS0.84H4.32"); // ss_AlFeSi-hydrogarnet
    SolidPhase ss_ettringite("ettringite ettringite30"); // ss_Ettrignite_ss
    SolidPhase ss_OH_SO4_AFm("C4AH13 monosulphate12"); // ss_Monosulfate_ss
    SolidPhase ss_CSHQ("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH"); // ss_C-S-H

    // Define chemical system by providing the database, aqueous phase, minerals phases (pure and solid solutions)
    ChemicalSystem system(db, aq_solution, pure_minerals, gaseous,
                          ss_C3AFS084H,
                          ss_ettringite,
                          ss_OH_SO4_AFm,
                          ss_CSHQ);

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

    // We define the materials for our equilibirum recipe
    // Cement clinker composition from XRF as given in Lothenbach et al., (2008) recalculated for 100g
    Material cement_clinker(system);
    cement_clinker.add("SiO2",  20.47, "g");
    cement_clinker.add("CaCO3", 0.00, "g");
    cement_clinker.add("CaO", 65.7, "g");
    cement_clinker.add("Al2O3", 4.9, "g");
    cement_clinker.add("Fe2O3", 3.2, "g");
    cement_clinker.add("K2O", 0.79, "g");
    cement_clinker.add("Na2O", 0.42, "g");
    cement_clinker.add("MgO", 1.8, "g");
    cement_clinker.add("SO3", 2.29, "g");
    cement_clinker.add("CO2", 0.26, "g");
    cement_clinker.add("O2", 0.15, "g");
    
    // Water
    Material water(system);
    water.add("H2O", 1000.0, "g");
    
    // We make a cement mix of 0.5 water/binder 
    Material cement_mix(system);
    cement_mix = cement_clinker(100.0, "g") + water(50.0, "g");
    
    // And we equilibrate our cement mix
    state = cement_mix.equilibrate(20.0, "celsius", 1.0, "bar");
    
    // Output the chemical state to a text-file
    state.output("state-cemdata18_1.txt");

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
    state.output("state-cemdata18_2.txt");

    // Output the chemical properties to a text-file
    ChemicalProps props(state);
    props.output("props.txt");

    // Output the aqueous properties to a text-file
    AqueousProps aprops(state);
    aprops.output("aprops.txt");

    return 0;
}
