// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

struct WaterHelmholtzState
{
    /// The specific Helmholtz free energy of water (in units of J/kg)
    real helmholtz = {};

    /// The first-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
    real helmholtzT = {};

    /// The first-order partial derivative of the specific Helmholtz free energy of water with respect to density
    real helmholtzD = {};

    /// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
    real helmholtzTT = {};

    /// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature and density
    real helmholtzTD = {};

    /// The second-order partial derivative of the specific Helmholtz free energy of water with respect to density
    real helmholtzDD = {};

    /// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
    real helmholtzTTT = {};

    /// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, temperature, and density
    real helmholtzTTD = {};

    /// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, density, and density
    real helmholtzTDD = {};

    /// The third-order partial derivative of the specific Helmholtz free energy of water with respect to density
    real helmholtzDDD = {};
};

} // namespace Reaktoro
