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
//   • Allan Leal (14 September 2021)
//   • William Smith (14 September 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database db = Database::fromFile(REAKTORO_EXAMPLES_DIR"/resources/db-molten-nacl.yaml");

    LiquidPhase liquid("NaCl(l)");

    GaseousPhase gases(speciate("Na Cl"));

    ChemicalSystem system(db, liquid, gases);

    ChemicalState state(system);
    state.temperature(1100, "K");
    state.pressure(100, "Pa");
    state.set("NaCl(l)", 1.0, "mol");

    equilibrate(state);

    ChemicalProps props(state);
    props.output("props.txt");

    std::cout << "Success! Check outputted file `props.txt`." << std::endl;

    return 0;
}
