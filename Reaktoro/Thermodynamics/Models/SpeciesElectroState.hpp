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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

struct SpeciesElectroState
{
    /// The effective electrostatic radius of the solute species at referente temperature 298.15 K and pressure 1 bar
    real reref = {};

    /// The effective electrostatic radius of the solute species
    real re = {};

    /// The Born coefficient of the solute species
    real w = {};

    /// The first-order partial derivative of the Born coefficient of the solute species with respect to temperature
    real wT = {};

    /// The first-order partial derivative of the Born coefficient of the solute species with respect to pressure
    real wP = {};

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to temperature
    real wTT = {};

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to temperature and pressure
    real wTP = {};

    /// The second-order partial derivative of the Born coefficient of the solute species with respect to pressure
    real wPP = {};
};

} // namespace Reaktoro
