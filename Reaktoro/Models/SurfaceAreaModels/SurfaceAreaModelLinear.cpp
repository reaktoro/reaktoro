// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "SurfaceAreaModelLinear.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

namespace Reaktoro {

auto SurfaceAreaModelLinearMolar(String const& phase, Param const& Abar) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        static thread_local auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).amount();
        return q * Abar;
    };
}

auto SurfaceAreaModelLinearSpecific(String const& phase, Param const& Abar) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        static thread_local auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).mass();
        return q * Abar;
    };
}

auto SurfaceAreaModelLinearVolumetric(String const& phase, Param const& Abar) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        static thread_local auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).volume();
        return q * Abar;
    };
}

auto SurfaceAreaModelLinear(String const& phase, real Abar, Chars unitAbar) -> SurfaceAreaModel
{
    if(units::convertible(unitAbar, "m2/mol"))
    {
        Abar = units::convert(Abar, unitAbar, "m2/mol");
        return SurfaceAreaModelLinearMolar(phase, Abar);
    }
    else if(units::convertible(unitAbar, "m2/kg"))
    {
        Abar = units::convert(Abar, unitAbar, "m2/kg");
        return SurfaceAreaModelLinearSpecific(phase, Abar);
    }
    else if(units::convertible(unitAbar, "m2/m3"))
    {
        Abar = units::convert(Abar, unitAbar, "m2/m3");
        return SurfaceAreaModelLinearVolumetric(phase, Abar);
    }
    else errorif(true, "Expecting a normalized surface area convertible to `m2/mol`, `m2/kg`, or `m2/m3` but got `", unitAbar, "` instead.")
}

} // namespace Reaktoro
