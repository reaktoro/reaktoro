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

#include "SurfaceAreaModelPower.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

namespace Reaktoro {

auto SurfaceAreaModelPowerMolar(String const& phase, real const& A0, real const& q0, real const& p) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        const auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).amount();
        return A0 * pow(q / q0, p);
    };
}

auto SurfaceAreaModelPowerSpecific(String const& phase, real const& A0, real const& q0, real const& p) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        const auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).mass();
        return A0 * pow(q / q0, p);
    };
}

auto SurfaceAreaModelPowerVolumetric(String const& phase, real const& A0, real const& q0, real const& p) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real
    {
        const auto iphase = props.system().phases().index(phase);
        const auto q = props.phaseProps(iphase).volume();
        return A0 * pow(q / q0, p);
    };
}

auto SurfaceAreaModelPower(String const& phase, real A0, Chars unitA0, real q0, Chars unitq0, real const& p) -> SurfaceAreaModel
{
    A0 = units::convert(A0, unitA0, "m2");

    if(units::convertible(unitq0, "mol"))
    {
        q0 = units::convert(q0, unitq0, "mol");
        return SurfaceAreaModelPowerMolar(phase, A0, q0, p);
    }
    else if(units::convertible(unitq0, "kg"))
    {
        q0 = units::convert(q0, unitq0, "kg");
        return SurfaceAreaModelPowerSpecific(phase, A0, q0, p);
    }
    else if(units::convertible(unitq0, "m3"))
    {
        q0 = units::convert(q0, unitq0, "m3");
        return SurfaceAreaModelPowerVolumetric(phase, A0, q0, p);
    }
    else errorif(true, "Expecting an initial phase quantity unit convertible to `mol`, `kg`, or `m3` but got `", unitq0, "` instead.")
}

} // namespace Reaktoro
