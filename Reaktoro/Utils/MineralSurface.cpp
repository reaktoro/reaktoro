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

#include "MineralSurface.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Models/SurfaceAreaModels/SurfaceAreaModelConstant.hpp>
#include <Reaktoro/Models/SurfaceAreaModels/SurfaceAreaModelLinear.hpp>
#include <Reaktoro/Models/SurfaceAreaModels/SurfaceAreaModelPower.hpp>

namespace Reaktoro {

MineralSurface::MineralSurface(String const& mineral)
: GeneralSurface(mineral) {}

MineralSurface::MineralSurface(String const& mineral, real A, Chars unitA)
: GeneralSurface(mineral)
{
    if(units::convertible(unitA, "m2"))
        setAreaModel(SurfaceAreaModelConstant(A, unitA));
    else try { setAreaModel(SurfaceAreaModelLinear(mineral, A, unitA)); }
    catch(...) {
        errorif(true, "Expecting surface area unit to be convertible to either `m2`, `m2/mol`, `m2/kg`, or `m2/m3`, but got `", unitA ,"` instead ");
    }
}

MineralSurface::MineralSurface(String const& mineral, real A0, Chars unitA0, real q0, Chars unitq0, real p)
: GeneralSurface(mineral, SurfaceAreaModelPower(mineral, A0, unitA0, q0, unitq0, p))
{}

} // namespace Reaktoro
