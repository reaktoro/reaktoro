// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "PhaseIdentification.hpp"

#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {
namespace PhaseIdentification {

auto volumeMethod(ThermoScalar Volume, ThermoScalar b) -> PhaseType
{
    if ((Volume.val / b.val) > 1.75)
        return PhaseType::Gas;

    return PhaseType::Liquid;
}

auto isothermalCompressibilityMethod(ThermoScalar Temperature, ThermoScalar Pressure, ChemicalScalar Z) -> PhaseType
{
    auto V = Z * universalGasConstant*T / P;
    auto dkdt = (1.0 / (V.val*V.val))*V.ddP*V.ddT;
    
    if (dkdt <= 0.0)
        return PhaseType::Gas;

    return PhaseType::Liquid;
}

auto 

}
}


