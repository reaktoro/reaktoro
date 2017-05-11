// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Gems gems("demos/resources/gems/CalciteBC-dat.lst");

    ChemicalSystem system = gems;

    ChemicalProperties properties(system);
    EquilibriumState state = gems.state(system);

    auto T = gems.temperature();
    auto P = gems.pressure();
    auto n = gems.speciesAmounts();
    properties.update(T, P, n);
    properties.update(T, P, n);
    properties.update(T, P, n);
    properties.update(T, P, n);
//
//    for(int i = 0; i < n.rows(); ++i) n[i] = std::max(n[i], 1e-10);

//    gems.set(T, P, n);
//    gems.set(T, P, n);
//    gems.set(T, P, n);
//    gems.set(T, P, n);

    std::cout << system << std::endl;
    std::cout << state << std::endl;

//    equilibrate(state);
}
