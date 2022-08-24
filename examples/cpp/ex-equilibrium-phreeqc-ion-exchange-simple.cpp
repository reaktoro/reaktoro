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
    AqueousPhase aqueous_phase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2");
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define an ion exchange phase
    IonExchangePhase exchange_phase("NaX KX CaX2 MgX2");
    exchange_phase.setActivityModel(ActivityModelIonExchangeGainesThomas());
    // Construct the chemical system
    ChemicalSystem system(db, aqueous_phase, exchange_phase);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in bar

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.setTemperature(T, "celsius");
    solutionstate.setPressure(P, "bar");
    solutionstate.setSpeciesMass("H2O"    , 1.00, "kg");
    solutionstate.setSpeciesAmount("Na+"  , 1.00, "mmol");
    solutionstate.setSpeciesAmount("Ca+2" , 1.00, "mmol");
    solutionstate.setSpeciesAmount("Mg+2" , 1.00, "mmol");
    solutionstate.setSpeciesAmount("K+"   , 1.00, "mmol");
    solutionstate.setSpeciesAmount("NaX"  , 1.00, "umol"); // set small to make sure we have plenty of water for available exchanger X-

    std::cout << "*******************************************" << std::endl;
    std::cout << "Before equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "b = " << solutionstate.elementAmounts().transpose() << std::endl;
    std::cout << "Z = " << solutionstate.charge() << std::endl;

    EquilibriumOptions opts;
    opts.optima.output.active = true;

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);

    auto res = solver.solve(solutionstate);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "b = " << solutionstate.elementAmounts().transpose() << std::endl;
    std::cout << "Z = " << solutionstate.charge() << std::endl;
    std::cout << "succeed = " << res.optima.succeeded << std::endl;
    std::cout << solutionstate << std::endl;

    AqueousProps aprops(solutionstate);
    std::cout << "I  = " << aprops.ionicStrength() << " mol/kgw" << std::endl;
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "pE = " << aprops.pE() << std::endl;

    //IonExchangeProps exprops(solutionstate);
    IonExchangeProps exprops(system);
    exprops.update(solutionstate);

    std::cout << exprops << std::endl;

    return 0;
}
