// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousSpecies;
struct SpeciesElectroState;
struct WaterThermoState;

/// A type used to describe the function g of the HKF model and its partial temperature and pressure derivatives
struct FunctionG
{
    /// The function g at temperature T and pressure P
    ThermoScalar g;

    /// The first-order partial derivative of function g with respect to temperature
    ThermoScalar gT;

    /// The first-order partial derivative of function g with respect to pressure
    ThermoScalar gP;

    /// The second-order partial derivative of function g with respect to temperature
    ThermoScalar gTT;

    /// The second-order partial derivative of function g with respect to temperature and pressure
    ThermoScalar gTP;

    /// The second-order partial derivative of function g with respect to pressure
    ThermoScalar gPP;
};

/// Calculate the function g of the HKF model.
auto functionG(ThermoScalar T, ThermoScalar P, const WaterThermoState& wts) -> FunctionG;

/// Calculate the electrostatic state of the aqueous species using the g-function state.
auto speciesElectroStateHKF(const FunctionG& g, const AqueousSpecies& species) -> SpeciesElectroState;

/// Calculate the electrostatic state of the aqueous species using the HKF model.
auto speciesElectroStateHKF(ThermoScalar T, ThermoScalar P, const AqueousSpecies& species) -> SpeciesElectroState;

} // namespace Reaktoro
