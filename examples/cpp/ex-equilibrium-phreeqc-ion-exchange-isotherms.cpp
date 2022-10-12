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
//   • Svetlana Kyas (23 November 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Core/Utils.hpp>

using namespace Reaktoro;

int main()
{
    // Initialize the Phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define an aqueous phase
    AqueousPhase aqueous_phase(speciate("H O C Ca Na Mg Cl K"));
    aqueous_phase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Define an ion exchange phase
    IonExchangePhase exchange_phase("KX CaX2");
    exchange_phase.setActivityModel(ActivityModelIonExchangeGainesThomas());

    // Construct the chemical system
    ChemicalSystem system(db, aqueous_phase, exchange_phase);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in bar

    // Initializing aqueous, chemical, and ion-exchange properties classes
    AqueousProps aprops(system);
    ChemicalProps chemprops(system);
    IonExchangeProps exprops(system);

    // To match PHREEQC example:
    // ------------------------------------------------------
    // SOLUTION 1
    // K 100.0
    // Ca 0.0
    //
    // EXCHANGE 1
    // KX 0.4 # 400 mmol KX # 400 mmol KX, total K = 500 mmol
    //
    // # Remove 500 mmol K+, replace with Ca2+...
    // REACTION; K -1 Ca 0.5; 0.500 in 20 steps
    // ------------------------------------------------------

    // Initialize auxiliary variables and K and Ca increments for the loop
    auto num_steps = 21;

    // Create an equilibrium solver
    EquilibriumSolver solver(system);

    // Output the header of the table
    std::cout << "   mols(K)   mols(Ca)          I   mols(K+) mols(Ca+2)   beta(KX) beta(CaX2)   mols(KX) mols(CaX2)" << std::endl;

    for(auto i = 0; i < num_steps; i++)
    {
        // The factor used below to regulate amounts of K+, KX, and Ca+2
        const auto factor = double(i)/(num_steps - 1);

        // Define initial chemical state
        ChemicalState solutionstate(system);
        solutionstate.temperature(T, "celsius");
        solutionstate.pressure(P, "atm");
        solutionstate.set("H2O"  , 1.00, "kg");
        solutionstate.set("K+"   , 0.10 * (1 - factor), "mol");
        solutionstate.set("KX"   , 0.40 * (1 - factor), "mol");
        solutionstate.set("Ca+2" , 0.25 * factor, "mol");

        auto res = solver.solve(solutionstate);

        errorif(res.failed(), "Calculation did not succeed!");

        aprops.update(solutionstate);
        chemprops.update(solutionstate);
        exprops.update(solutionstate);

        std::cout.precision(4);
        std::cout << std::scientific
            << chemprops.elementAmount("K") << " "
            << chemprops.elementAmount("Ca") << " "
            << aprops.ionicStrength() << " "
            << chemprops.speciesAmount("K+") << " "
            << chemprops.speciesAmount("Ca+2") << " "
            << exprops.speciesEquivalentFraction("KX") << " "
            << exprops.speciesEquivalentFraction("CaX2") << " "
            << exprops.speciesAmount("KX") << " "
            << exprops.speciesAmount("CaX2") << std::endl;
    }

    return 0;
}
