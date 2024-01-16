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

#include "ReactionThermoProps.hpp"

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

auto operator<<(std::ostream& out, const ReactionThermoProps& props) -> std::ostream&
{
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(props.T), "K" });
    table.add_row({ "Pressure", strfix(props.P*1e-5), "bar" });
    table.add_row({ "Equilibrium Constant (log base 10)", strfix(props.lgK), "-" });
    table.add_row({ "Delta Standard Gibbs Energy", strfix(props.dG0), "J/mol" });
    table.add_row({ "Delta Standard Enthalpy", strfix(props.dH0), "J/mol" });
    table.add_row({ "Delta Standard Volume", strsci(props.dV0), "m3/mol" });
    table.add_row({ "Delta Standard Volume (Temperature Derivative)", strsci(props.dVT0), "m3/(mol*K)" });
    table.add_row({ "Delta Standard Volume (Pressure Derivative)", strsci(props.dVP0), "m3/(mol*Pa)" });
    table.add_row({ "Delta Standard Isobaric Heat Capacity", strfix(props.dCp0), "J/(mol*K)" });
    table.add_row({ "Delta Standard Isochoric Heat Capacity", strfix(props.dCv0), "J/(mol*K)" });
    table.add_row({ "Delta Standard Internal Energy", strfix(props.dU0), "J/mol" });
    table.add_row({ "Delta Standard Entropy", strfix(props.dS0), "J/(mol*K)" });
    table.add_row({ "Delta Standard Helmholtz Energy", strfix(props.dA0), "J/mol" });

    auto i = 0;
    for(auto& row : table)
    {
        if(i >= 2)  // apply from the third row
            table[i]
                .format()
                .border_top("")
                .column_separator("")
                .corner_top_left("")
                .corner_top_right("");
        i += 1;
    }

    table.row(0).format().font_style({FontStyle::bold});  // Bold face for header
    table.column(1).format().font_align(FontAlign::right); // Value column with right alignment
    table.column(2).format().font_align(FontAlign::right); // Unit column with right alignment

    auto old_locale = std::locale::global(std::locale("C")); // This locale logic is needed to avoid UnicodeDecodeError: 'utf-8' codec can't decode byte 0xa0 in position ...
    out << table;
    std::locale::global(old_locale);

    return out;
}

} // namespace Reaktoro
