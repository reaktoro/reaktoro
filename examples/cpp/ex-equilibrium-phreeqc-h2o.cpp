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
    PhreeqcDatabase db("phreeqc.dat");

    Phases phases(db,  AqueousSolution("H2O H+ OH- O2 H2"));

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

    const auto H2O = system.species(0);
    const auto Hp = system.species(1);
    const auto OHm = system.species(2);
    const auto O2 = system.species(3);
    const auto H2 = system.species(4);

    const auto u0_H2O = H2O.props(T, P).G0;
    const auto u0_Hp  = Hp.props(T, P).G0;
    const auto u0_OHm = OHm.props(T, P).G0;
    const auto u0_O2  = O2.props(T, P).G0;
    const auto u0_H2  = H2.props(T, P).G0;

    std::cout << "u0_H2O = " << u0_H2O << std::endl;
    std::cout << "u0_Hp  = " << u0_Hp << std::endl;
    std::cout << "u0_OHm = " << u0_OHm << std::endl;
    std::cout << "u0_O2  = " << u0_O2 << std::endl;
    std::cout << "u0_H2  = " << u0_H2 << std::endl;



    // ChemicalState state(system);
    // // state.setTemperature(60.0, "celsius");
    // // state.setPressure(30.0, "bar");
    // state.setTemperature(25.0, "celsius");
    // state.setPressure(1.0, "bar");
    // state.setSpeciesMass("H2O", 1.0, "kg");

    // EquilibriumOptions options;
    // options.optima.output.active = true;

    // EquilibriumSolver solver(system);
    // solver.setOptions(options);

    // solver.solve(state);

    // std::cout << state.speciesAmounts() << std::endl;

    return 0;
}
