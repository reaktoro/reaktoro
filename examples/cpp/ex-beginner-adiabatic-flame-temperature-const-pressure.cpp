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
//   • Allan Leal (19 July 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <iostream>
using namespace std;

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    String reactantName = "CH4";
    String oxidizerName = "Air";

    double reactantAmount = 1.0; // in mol
    double oxidizerAmount = 2.0; // in mol

    cout << "LOADING NASA THERMO.INP DATABASE...";
    NasaDatabase db("nasa-cea");

    Species reactant = db.species().getWithName(reactantName);
    Species oxidizer = db.species().getWithName(oxidizerName);

    Strings symbols = unique(concatenate(reactant.elements().symbols(), oxidizer.elements().symbols()));

    GaseousPhase gases(speciate(symbols));
    CondensedPhases condensed(speciate(symbols));

    cout << "CREATING CHEMICAL SYSTEM..." << endl;

    ChemicalSystem system(db, gases, condensed);

    cout << "CHEMICAL SYSTEM CREATED WITH " << system.elements().size() << " ELEMENTS AND " << system.species().size() << " SPECIES" << endl;

    ChemicalState state0(system);
    state0.temperature(25.0, "celsius");
    state0.pressure(1.0, "atm");
    state0.setSpeciesAmounts(1e-16);
    state0.set(reactantName, reactantAmount, "mol");
    state0.set(oxidizerName, oxidizerAmount, "mol");

    ChemicalProps props0(state0);

    EquilibriumSpecs specs(system);
    specs.pressure();
    specs.enthalpy();

    EquilibriumConditions conditions(specs);
    conditions.pressure(props0.pressure());
    conditions.enthalpy(props0.enthalpy());
    conditions.setLowerBoundTemperature(298.15, "celsius");
    conditions.setUpperBoundTemperature(4000.0, "celsius");

    ChemicalState state(state0);

    cout << "PERFORMING THE CHEMICAL EQUILIBRIUM CALCULATION..." << endl;
    EquilibriumSolver solver(specs);

    auto result = solver.solve(state, conditions);

    if(result.optima.succeeded)
    {
        cout << state << endl;
        cout << "COMPUTED ADIABATIC FLAME TEMPERATURE: " << state.temperature() << " K" << endl;
        cout << "ITERATIONS: " << result.optima.iterations << endl;
    }
    else cout << "ERROR!" << endl;

    return 0;
}
