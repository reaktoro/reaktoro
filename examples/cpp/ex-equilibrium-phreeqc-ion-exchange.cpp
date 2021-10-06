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
//   • Svetlana Kyas (30 August 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Core/Utils.hpp>

using namespace Reaktoro;

const auto T = 25.0; // temperature in celsius
const auto P = 1.0;  // pressure in bar

int main()
{
    // Initialize the Phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define an aqueous phase
    AqueousPhase aqueous_phase(speciate("H O C Ca Na Mg Cl"));
    aqueous_phase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Fetch species for ion-exchange modeling
    SpeciesList exchange_species = db.species().withAggregateState(AggregateState::IonExchange);

    // Define an ion exchange phase
    IonExchangePhase exchange_phase(detail::extractNames(exchange_species));
    //IonExchangePhase exchange_phase(speciate("X Ca Na Mg"));
    //IonExchangePhase exchange_phase(speciate("X Ca Na Mg"), exclude("organic")); //
    exchange_phase.setActivityModel(ActivityModelIonExchangeGainesThomas());

    // Construct the chemical system
    ChemicalSystem system(db, aqueous_phase, exchange_phase);

    // Specify conditions to be satisfied at the chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    //specs.charge();
    specs.pH();
    //specs.openTo("Cl-");

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");
    //conditions.charge(0.0);
    conditions.pH(7.0);

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.setSpeciesMass("H2O"    , 1.00, "kg");
//    solutionstate.setSpeciesAmount("Na+"  , 1.10, "mol");
//    solutionstate.setSpeciesAmount("Mg+2" , 0.48, "mol");
//    solutionstate.setSpeciesAmount("Ca+2" , 1.90, "mol");
    solutionstate.setSpeciesAmount("Na+"  , 1.10, "mmol");
    solutionstate.setSpeciesAmount("Mg+2" , 0.48, "mmol");
    solutionstate.setSpeciesAmount("Ca+2" , 1.90, "mmol");
    solutionstate.setSpeciesAmount("X-"   , 0.06, "mol");

    // Output the chemical state to a text file
    std::cout << solutionstate << std::endl;

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);
    solver.solve(solutionstate, conditions);

    // Output the chemical state to a text file
    solutionstate.output("state.txt");
    std::cout << solutionstate << std::endl;

    AqueousProps aprops(solutionstate);
    std::cout << "I  = " << aprops.ionicStrength() << " mol/kgw" << std::endl;
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "pE = " << aprops.pE() << std::endl;

    return 0;
}
