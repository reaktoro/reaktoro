// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Core/Utils.hpp>

using namespace Reaktoro;

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

    // Define an ion exchange phase
    IonExchangePhase exchange_phase("NaX CaX2 MgX2");
    exchange_phase.setActivityModel(ActivityModelIonExchangeGainesThomas());

    // Construct the chemical system
    ChemicalSystem system(db, aqueous_phase, exchange_phase);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in atm

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.setTemperature(T, "celsius");
    solutionstate.setPressure(P, "atm");
    solutionstate.setSpeciesMass("H2O"    , 1.00, "kg");
    // Scale solution recipe to match the values of the PHREEQC examples
    solutionstate.setSpeciesAmount("Na+"  , 1.10, "mol");
    solutionstate.setSpeciesAmount("Mg+2" , 0.48, "mol");
    solutionstate.setSpeciesAmount("Ca+2" , 1.90, "mol");
    // Set the number of exchange assuming that it is completely occupied by Na
    solutionstate.setSpeciesAmount("NaX"  , 0.06, "umol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(solutionstate);

    // Output the chemical state to a text file
    solutionstate.output("state.txt");
    std::cout << solutionstate << std::endl;

    AqueousProps aprops(solutionstate);
    std::cout << "I  = " << aprops.ionicStrength() << " mol/kgw" << std::endl;
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "pE = " << aprops.pE() << std::endl;

    ChemicalProps chemprops(solutionstate);
    IonExchangeProps exprops(chemprops);
    std::cout << exprops << std::endl;

    return 0;
}
