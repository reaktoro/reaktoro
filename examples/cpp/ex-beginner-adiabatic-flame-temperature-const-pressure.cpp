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

const auto xN2  = 0.78084;
const auto xO2  = 0.209476;
const auto xAr  = 0.9365;
const auto xCO2 = 0.0319;

int main()
{
    // String reactantName;
    // String oxidizerName;
    // double reactantAmount = 1.0; // in mol
    // double oxidizerAmount = 1.0; // in mol

    // cout << "ENTER NAME OF REACTANT SPECIES:" << endl; cin >> reactantName;
    // cout << "ENTER NAME OF OXIDIZER SPECIES:" << endl; cin >> oxidizerName;

    // cout << "ENTER AMOUNT OF REACTANT (IN MOL):" << endl; cin >> reactantAmount;
    // cout << "ENTER AMOUNT OF OXIDIZER (IN MOL):" << endl; cin >> oxidizerAmount;

    String reactantName = "Mg(cr)";
    // String oxidizerName = "Air";
    String oxidizerName = "O2";
    double reactantAmount = 1.0; // in mol
    double oxidizerAmount = 0.5; // in mol

    cout << "LOADING NASA THERMO.INP DATABASE...";
    NasaDatabase db("cea");

    Species reactant = db.species().getWithName(reactantName);
    Species oxidizer = db.species().getWithName(oxidizerName);

    Strings symbols = unique(concatenate(reactant.elements().symbols(), oxidizer.elements().symbols()));
    symbols.push_back("N");
    symbols.push_back("Ar");
    symbols.push_back("C");

    GaseousPhase gases(speciate(symbols));
    // GaseousPhase gases("e- Ar Ar+ C C+ C- CN CN+ CN- CNN CO CO+ CO2 CO2+ C2 C2+ C2- CCN CNC OCCN C2N2 C2O C3 CNCOCN C3O2 C4 C4N2 C5 Mg Mg+ MgN MgO Mg2 N N+ N- NCO NO NO+ NO2 NO2- NO3 NO3- N2 N2+ N2- NCN N2O N2O+ N2O3 N2O4 N2O5 N3 O O+ O- O2 O2+ O2- O3");
    // GaseousPhase gases("Ar C CN CNN CO CO2 C2 CCN CNC OCCN C2N2 C2O C3 CNCOCN C3O2 C4 C4N2 C5 Mg MgN MgO Mg2 N NCO NO NO2 NO3 N2 NCN N2O N2O3 N2O4 N2O5 N3 O O2 O3");
    // GaseousPhase gases("Ar CO CO2 Mg MgO Mg2 N2 O2");
    // GaseousPhase gases("C CN CNN CO CO2 C2 CCN CNC OCCN C2N2 C2O C3 CNCOCN C3O2 C4 C4N2 C5 Mg MgN MgO Mg2 N NCO NO NO2 NO3 N2 NCN N2O N2O3 N2O4 N2O5 N3 O O2 O3 Air");
    // GaseousPhase gases("Mg MgO Mg2 O O2 O3");

    Phases phases(db);
    phases.add(gases);

    // phases.add(CondensedPhase("Mg(cr)"));
    // phases.add(CondensedPhase("MgO(cr)"));
    // phases.add(CondensedPhase("Mg(L)"));

    // phases.add(CondensedPhase("Mg(cr)"));
    // phases.add(CondensedPhase("Mg(L)"));
    // phases.add(CondensedPhase("MgO(cr)"));
    // phases.add(CondensedPhase("MgO(L)"));
    // phases.add(CondensedPhase("MgCO3(cr)"));
    // phases.add(CondensedPhase("MgCO3(L)"));

    // phases.add(CondensedPhase("C(gr)"));
    // phases.add(CondensedPhase("Mg(cr)"));
    // phases.add(CondensedPhase("Mg(L)"));
    // phases.add(CondensedPhase("MgCO3(cr)"));
    // phases.add(CondensedPhase("MgCO3(L)"));
    // phases.add(CondensedPhase("MgO(cr)"));
    // phases.add(CondensedPhase("MgO(L)"));
    // phases.add(CondensedPhase("Mg3N2(cr)"));
    // phases.add(CondensedPhase("C2N2(L)"));
    // phases.add(CondensedPhase("N2(L)"));
    // phases.add(CondensedPhase("N2O4(L)"));
    // phases.add(CondensedPhase("O2(L)"));
    // phases.add(CondensedPhase("O3(L)"));

    // phases.add(CondensedPhase("C(gr)"));
    // phases.add(CondensedPhase("Mg(cr)"));
    // phases.add(CondensedPhase("Mg(L)"));
    // phases.add(CondensedPhase("MgCO3(cr)"));
    // phases.add(CondensedPhase("MgCO3(L)"));
    // phases.add(CondensedPhase("MgO(cr)"));
    // phases.add(CondensedPhase("MgO(L)"));
    // phases.add(CondensedPhase("Mg3N2(cr)"));
    // phases.add(CondensedPhase("C2N2(L)"));
    // phases.add(CondensedPhase("N2(L)"));
    // phases.add(CondensedPhase("N2O4(L)"));
    // phases.add(CondensedPhase("O2(L)"));
    // phases.add(CondensedPhase("O3(L)"));


    // phases.add(CondensedPhase("H2O(L)"));
    // phases.add(CondensedPhase("H2O(cr)"));
    // phases.add(CondensedPhase("H2O2(L)"));
    // phases.add(CondensedPhase("H2(L)"));

    // if(reactant.aggregateState() == AggregateState::CondensedPhase)
    //     phases.add(CondensedPhase(reactantName));

    // phases.add(CondensedPhase("MgO(cr)"));

    // phases.add(CondensedPhase("H2O(cr)"));
    // phases.add(CondensedPhase("H2O(L)"));
    // phases.add(CondensedPhase("H2(L)"));
    // phases.add(CondensedPhase("H2O2(L)"));
    // phases.add(CondensedPhase("O2(L)"));
    // phases.add(CondensedPhase("O3(L)"));

    CondensedPhases condensed(speciate(symbols));
    phases.add(condensed);

    cout << "CREATING CHEMICAL SYSTEM..." << endl;

    ChemicalSystem system(phases);

    cout << "CHEMICAL SYSTEM CREATED WITH " << system.elements().size() << " ELEMENTS AND " << system.species().size() << " SPECIES" << endl;

    ChemicalState state0(system);
    state0.temperature(25.0, "celsius");
    state0.pressure(1.0, "atm");
    state0.setSpeciesAmounts(1e-16);
    state0.set(reactantName, reactantAmount, "mol");
    state0.set(oxidizerName, oxidizerAmount, "mol");
    state0.set("N2", oxidizerAmount/xO2 * xN2, "mol");
    state0.set("CO2", oxidizerAmount/xO2 * xCO2, "mol");
    state0.set("Ar", oxidizerAmount/xO2 * xAr, "mol");

    ChemicalProps props0(state0);

    EquilibriumSpecs specs(system);
    specs.pressure();
    specs.enthalpy();

    EquilibriumConditions conditions(specs);
    conditions.pressure(props0.pressure());
    conditions.enthalpy(props0.enthalpy());
    // conditions.setLowerBoundTemperature(1000.0, "celsius");
    // conditions.setUpperBoundTemperature(4000.0, "celsius");
    conditions.setLowerBoundTemperature(0.0, "celsius");
    conditions.setUpperBoundTemperature(4000.0, "celsius");

    ChemicalState state(state0);

    EquilibriumOptions options;
    options.optima.output.active = true;
    options.optima.maxiterations = 100;
    options.logarithm_barrier_factor = 1.0e4; // IMPORTANT: Because 999'999'999'999 is assigned to G0 when temperature is out-of-range, the log barrier factor needs to be higher. If epsilon is 1e-16, a factor of about 1e4 is needed to make tau = 1e-12 instead of tau = 1e-16. Otherwise, very large Newton steps are generated and mass conservation strongly deviated.

    cout << "PERFORMING THE CHEMICAL EQUILIBRIUM CALCULATION..." << endl;
    EquilibriumSolver solver(specs);
    solver.setOptions(options);


    auto result = solver.solve(state, conditions);


    options.logarithm_barrier_factor = 1.0;
    solver.setOptions(options);
    result = solver.solve(state, conditions);

    if(result.optima.succeeded)
    {
        cout << state << endl;
        cout << "COMPUTED ADIABATIC FLAME TEMPERATURE: " << state.temperature() << " K" << endl;
        cout << "ITERATIONS: " << result.optima.iterations << endl;
    }
    else cout << "ERROR!" << endl;


    return 0;
}
