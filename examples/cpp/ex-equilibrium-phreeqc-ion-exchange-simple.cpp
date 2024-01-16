// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
//   • Allan Leal (28 August 2023)
//     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>

using namespace Reaktoro;

int main()
{
    // Initialize the Phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define an aqueous phase
    AqueousPhase aqueousphase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2");
    aqueousphase.set(ActivityModelPhreeqc(db));

    // Define an ion exchange phase
    IonExchangePhase exchangephase("NaX KX CaX2 MgX2");
    exchangephase.set(ActivityModelIonExchangeGainesThomas());

    // Construct the chemical system
    ChemicalSystem system(db, aqueousphase, exchangephase);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in bar

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(T, "celsius");
    solutionstate.pressure(P, "bar");
    solutionstate.set("H2O"    , 1.00, "kg");
    solutionstate.set("Na+"  , 1.00, "mmol");
    solutionstate.set("Ca+2" , 1.00, "mmol");
    solutionstate.set("Mg+2" , 1.00, "mmol");
    solutionstate.set("K+"   , 1.00, "mmol");
    solutionstate.set("NaX"  , 1.00, "umol"); // set small to make sure we have plenty of water for available exchanger X-

    std::cout << "*******************************************" << std::endl;
    std::cout << "Before equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "b = " << solutionstate.elementAmounts().transpose() << std::endl;
    std::cout << "Z = " << solutionstate.charge() << std::endl;

    EquilibriumOptions opts;
    opts.optima.output.active = true;

    // Create an equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);

    auto res = solver.solve(solutionstate);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "b = " << solutionstate.elementAmounts().transpose() << std::endl;
    std::cout << "Z = " << solutionstate.charge() << std::endl;
    std::cout << "succeed = " << res.succeeded() << std::endl;
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
