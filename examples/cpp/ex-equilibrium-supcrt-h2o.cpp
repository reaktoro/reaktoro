// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    SupcrtDatabase db("supcrt98.xml");

    AqueousSolution aqueousphase(speciate("H O C Na Cl"));
    // AqueousSolution aqueousphase("H2O H+ OH- HCO3- CO2");
    aqueousphase.setActivityModel(ActivityModelHKF());

    GaseousSolution gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    // Phases phases(db, AqueousSolution("H2O H+ OH- O2 H2"));
    // Phases phases(db,  AqueousSolution("H2O H+ OH-"));
    Phases phases(db, aqueousphase, gaseousphase);

    ChemicalSystem system(db, phases);

    const auto A = system.formulaMatrix();

    for(auto s : system.species())
    {
        std::cout << s.name() << std::endl;

        for(auto [element, coeff] : s.elements())
            std::cout << element.symbol() << "*" << coeff << std::endl;
    }


    for(auto e : system.elements())
        std::cout << e.symbol() << std::endl;

    std::cout << "A = \n" << A << std::endl;

    double T = 25.0 + 273.15;
    double P = 1.0 * 1e5;

    // const auto H2O = system.species(0);
    // const auto Hp = system.species(1);
    // const auto OHm = system.species(2);
    // const auto O2 = system.species(3);
    // const auto H2 = system.species(4);

    // const auto u0_H2O = H2O.props(T, P).G0;
    // const auto u0_Hp  = Hp.props(T, P).G0;
    // const auto u0_OHm = OHm.props(T, P).G0;
    // const auto u0_O2  = O2.props(T, P).G0;
    // const auto u0_H2  = H2.props(T, P).G0;

    // std::cout << "u0_H2O = " << u0_H2O << std::endl;
    // std::cout << "u0_Hp  = " << u0_Hp << std::endl;
    // std::cout << "u0_OHm = " << u0_OHm << std::endl;
    // std::cout << "u0_O2  = " << u0_O2 << std::endl;
    // std::cout << "u0_H2  = " << u0_H2 << std::endl;

    const auto H2O = system.species(0);
    const auto Hp = system.species(1);
    const auto OHm = system.species(2);
    const auto HCO3m = system.species(3);
    const auto CO2 = system.species(4);
    const auto CO2g = system.species(5);

    const auto G0_H2O = H2O.props(T, P).G0;
    const auto G0_Hp = Hp.props(T, P).G0;
    const auto G0_OHm = OHm.props(T, P).G0;
    const auto G0_HCO3m = HCO3m.props(T, P).G0;
    const auto G0_CO2 = CO2.props(T, P).G0;
    const auto G0_CO2g = CO2g.props(T, P).G0;

    std::cout << "G0_H2O   = " << G0_H2O << std::endl;
    std::cout << "G0_Hp    = " << G0_Hp << std::endl;
    std::cout << "G0_OHm   = " << G0_OHm << std::endl;
    std::cout << "G0_HCO3m = " << G0_HCO3m << std::endl;
    std::cout << "G0_CO2   = " << G0_CO2 << std::endl;
    std::cout << "G0_CO2g  = " << G0_CO2g << std::endl;

    const auto R = universalGasConstant;

    std::cout << "lgK_CO2g  = " << -(G0_CO2g - G0_CO2)/(R*T*ln10) << std::endl;




    ChemicalState state(system);
    // state.setTemperature(60.0, "celsius");
    // state.setPressure(30.0, "bar");
    // state.setTemperature(60.0, "celsius");
    // state.setPressure(100.0, "bar");
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    // state.setSpeciesMass("H2O", 1.0, "kg");
    state.setSpeciesAmount("H2O(l)", 55.0, "mol");
    state.setSpeciesAmount("CO2(g)", 50.0, "mol");
    state.setSpeciesAmount("Na+", 1.0, "mol");
    state.setSpeciesAmount("Cl-", 1.0, "mol");
    // state.setSpeciesAmount("O2", 1.0e-14, "mol"); // FIXME: Without this, the example currently fails.
    // state.setSpeciesAmount("O2", 1.0e-16, "mol"); // FIXME: Without this, the example currently fails.
    // state.setSpeciesAmount("O2", 1.0e-35, "mol"); // FIXME: Without this, the example currently fails.
    // state.setSpeciesAmount("O2", 0.999999999e-39, "mol"); // FIXME: Without this, the example currently fails.

    EquilibriumOptions options;
    options.optima.output.active = true;
    options.optima.max_iterations = 50;
    // options.optima.kkt.method = Optima::SaddlePointMethod::Fullspace;
    // options.optima.kkt.method = Optima::SaddlePointMethod::Nullspace;
    options.optima.kkt.method = Optima::SaddlePointMethod::Rangespace;
    options.optima.linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor = 100000000000.0;
    options.optima.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = 100000000000.0;

    EquilibriumSolver solver(system);
    solver.setOptions(options);

    solver.solve(state);

    const auto n = state.speciesAmounts();

    for(auto i = 0; i < n.size(); ++i)
    {
        std::cout << std::setw(20) << system.species(i).name();
        std::cout << std::setw(20) << n[i];
        std::cout << std::endl;
    }



    return 0;
}
