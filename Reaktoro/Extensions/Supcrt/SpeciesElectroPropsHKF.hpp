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

#pragma once

// Reaktoro includes
#include <Reaktoro/Extensions/Supcrt/SupcrtParams.hpp>

namespace Reaktoro {

// Forward declarations
struct SpeciesElectroProps;
struct WaterThermoProps;

/// A type used to describe the function g of the HKF model and its partial temperature and pressure derivatives
struct FunctionG
{
    /// The function g at temperature T and pressure P
    real g = {};

    /// The first-order partial derivative of function g with respect to temperature
    real gT = {};

    /// The first-order partial derivative of function g with respect to pressure
    real gP = {};

    /// The second-order partial derivative of function g with respect to temperature
    real gTT = {};

    /// The second-order partial derivative of function g with respect to temperature and pressure
    real gTP = {};

    /// The second-order partial derivative of function g with respect to pressure
    real gPP = {};
};

/// Calculate the function g of the HKF model.
auto functionG(real T, real P, const WaterThermoProps& wtp) -> FunctionG;

/// Calculate the electrostatic state of the aqueous species using the g-function state.
auto speciesElectroPropsHKF(const FunctionG& g, const SupcrtParamsAqueousSoluteHKF& params) -> SpeciesElectroProps;

/// Calculate the electrostatic state of the aqueous species using the HKF model.
auto speciesElectroPropsHKF(real T, real P, const SupcrtParamsAqueousSoluteHKF& params) -> SpeciesElectroProps;

} // namespace Reaktoro
