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

#include "SpeciesThermoProps.hpp"

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

SpeciesThermoProps::SpeciesThermoProps(const real& T, const real& P, const StandardThermoProps& sprops)
: T(T),
  P(P),
  G0(sprops.G0),
  H0(sprops.H0),
  V0(sprops.V0),
  VT0(sprops.VT0),
  VP0(sprops.VP0),
  Cp0(sprops.Cp0),
  Cv0(VP0 == 0.0 ? Cp0 : Cp0 + T*VT0*VT0/VP0),
  U0(H0 - P*V0),
  S0((H0 - G0)/T),
  A0(U0 - T*S0)
{}

auto operator<<(std::ostream& out, const SpeciesThermoProps& props) -> std::ostream&
{
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", str(props.T), "K" });
    table.add_row({ "Pressure", str(props.P), "Pa" });
    table.add_row({ "Standard Gibbs Energy", str(props.G0), "J/mol" });
    table.add_row({ "Standard Enthalpy", str(props.H0), "J/mol" });
    table.add_row({ "Standard Volume", str(props.V0), "m3/mol" });
    table.add_row({ "Standard Volume (Temperature Derivative)", str(props.VT0), "m3/(mol*K)" });
    table.add_row({ "Standard Volume (Pressure Derivative)", str(props.VP0), "m3/(mol*Pa)" });
    table.add_row({ "Standard Isobaric Heat Capacity", str(props.Cp0), "J/(mol*K)" });
    table.add_row({ "Standard Isochoric Heat Capacity", str(props.Cv0), "J/(mol*K)" });
    table.add_row({ "Standard Internal Energy", str(props.U0), "J/mol" });
    table.add_row({ "Standard Entropy", str(props.S0), "J/(mol*K)" });
    table.add_row({ "Standard Helmholtz Energy", str(props.A0), "J/mol" });

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

    out << table;
    return out;
}

auto operator+(const SpeciesThermoProps& l, const SpeciesThermoProps& r) -> SpeciesThermoProps
{
    errorif(l.T != r.T, "Expecting identical temperatures when adding two SpeciesThermoProps objects.");
    errorif(l.P != r.P, "Expecting identical pressures when adding two SpeciesThermoProps objects.");

    return SpeciesThermoProps{
        l.T,
        l.P,
        l.G0  + r.G0,
        l.H0  + r.H0,
        l.V0  + r.V0,
        l.VT0 + r.VT0,
        l.VP0 + r.VP0,
        l.Cp0 + r.Cp0,
        l.Cv0 + r.Cv0,
        l.U0  + r.U0,
        l.S0  + r.S0,
        l.A0  + r.A0,
    };
}

auto operator-(const SpeciesThermoProps& l, const SpeciesThermoProps& r) -> SpeciesThermoProps
{
    errorif(l.T != r.T, "Expecting identical temperatures when subtracting two SpeciesThermoProps objects.");
    errorif(l.P != r.P, "Expecting identical pressures when subtracting two SpeciesThermoProps objects.");

    return SpeciesThermoProps{
        l.T,
        l.P,
        l.G0  - r.G0,
        l.H0  - r.H0,
        l.V0  - r.V0,
        l.VT0 - r.VT0,
        l.VP0 - r.VP0,
        l.Cp0 - r.Cp0,
        l.Cv0 - r.Cv0,
        l.U0  - r.U0,
        l.S0  - r.S0,
        l.A0  - r.A0,
    };
}

auto operator*(double l, const SpeciesThermoProps& r) -> SpeciesThermoProps
{
    return SpeciesThermoProps{
        r.T,
        r.P,
        r.G0  * l,
        r.H0  * l,
        r.V0  * l,
        r.VT0 * l,
        r.VP0 * l,
        r.Cp0 * l,
        r.Cv0 * l,
        r.U0  * l,
        r.S0  * l,
        r.A0  * l,
    };
}

auto operator*(const SpeciesThermoProps& l, double r) -> SpeciesThermoProps
{
    return r * l;
}

auto operator/(const SpeciesThermoProps& l, double r) -> SpeciesThermoProps
{
    return SpeciesThermoProps{
        l.T,
        l.P,
        l.G0  / r,
        l.H0  / r,
        l.V0  / r,
        l.VT0 / r,
        l.VP0 / r,
        l.Cp0 / r,
        l.Cv0 / r,
        l.U0  / r,
        l.S0  / r,
        l.A0  / r,
    };
}

} // namespace Reaktoro
