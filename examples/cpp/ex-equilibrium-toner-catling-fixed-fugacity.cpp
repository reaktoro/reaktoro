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
//   • Svetlana Kyas (14 July 2021)
//
// and since revised by:
//   • Allan Leal (16 July 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Initialize a thermodynamic database
    PhreeqcDatabase db = PhreeqcDatabase::fromFile(REAKTORO_EXAMPLES_DIR"/resources/phreeqc-toner-catling.dat");

    // Define aqueous phase
    AqueousPhase aqueousphase(speciate("H O C Na Cl"));
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Define gaseous phase
    MineralPhases minerals("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O");

    // Construct the chemical system
    ChemicalSystem system(db, aqueousphase, minerals);
    ChemicalProps props(system);
    AqueousProps aprops(system);

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.fugacity("CO2");

    EquilibriumConditions conditions(specs);

    auto T = 50.0; // temperature in celsius
    auto P = 1.0;  // pressure in bar

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);

    EquilibriumOptions opts;
    opts.epsilon = 1e-13;
    opts.optima.output.active = true;
    solver.setOptions(opts);

    EquilibriumResult res;

    // Define initial equilibrium state
    ChemicalState state0(system);
    state0.set("H2O"        ,   1.00, "kg");
    state0.set("Nahcolite"  ,  10.00, "mol");    // NaHCO3
    state0.set("Natron"     ,   0.00, "mol");    // Na2CO3:10H2O
    state0.set("Trona"      ,   0.00, "mol");    // Na3H(CO3)2:2H2O
    state0.set("Na2CO3:H2O" ,   0.00, "mol");
    state0.set("Na2CO3:7H2O",   0.00, "mol");
    state0.set("CO2"        , 100.00, "mol");

    // Initialize auxiliary variables and K and Ca increments for the loop
    auto num_ppco2s = 71; auto initialize = true;  auto d_ppCO2 = 0.1; //ppCO2 = -2.3 and -2.4, does not converge for T = 50C (if state is initialized in the loop)

    // Initial partial pressure of CO2
    auto ppCO2 = 2.0;

    // Output the header of the table
    std::cout << "     ppCO2   succeeded         pH  mols(CO3-2)  mols(HCO3-)" << std::endl;


    for(int i = 0; i < num_ppco2s; i++) {

        if((ppCO2 > -2.6) && (ppCO2 <= -2.2))
        //if(true)
        {
            ChemicalState state = state0;

            conditions.temperature(T, "celsius");
            conditions.pressure(P, "atm");
            conditions.fugacity("CO2", pow(10, ppCO2), "atm");

            if(initialize)
            {
                state = ChemicalState(system);
                state.set("H2O"        ,   1.00, "kg");
                state.set("Nahcolite"  ,  10.00, "mol");    // NaHCO3
                state.set("Natron"     ,   0.00, "mol");    // Na2CO3:10H2O
                state.set("Trona"      ,   0.00, "mol");    // Na3H(CO3)2:2H2O
                state.set("Na2CO3:H2O" ,   0.00, "mol");
                state.set("Na2CO3:7H2O",   0.00, "mol");
                state.set("CO2"        , 100.00, "mol");
            }

            res = solver.solve(state, conditions);

            aprops.update(state);
            props.update(state);

            std::cout.precision(4);
            std::cout << std::scientific
                      << ppCO2 << "           "
                      << res.optima.succeeded << " "
                      << aprops.pH() << "   "
                      << state.speciesAmount("CO3-2") << "   "
                      << state.speciesAmount("HCO3-") << std::endl;
        }
        ppCO2 -= d_ppCO2;
    }

    return 0;
}
