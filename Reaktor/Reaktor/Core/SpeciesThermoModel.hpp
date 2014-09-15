// Reaktor is a C++ library for computational reaction modelling.
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

// C++ includes
#include <functional>

// Reaktor includes
#include <Reaktor/Common/ThermoScalar.hpp>

namespace Reaktor {

/// Describe the thermodynamic model of a species
struct SpeciesThermoModel
{
    /// The standard molar volume function of the species (in units of m3/mol).
    std::function<ThermoScalar(double, double)> V;

    /// The standard molar entropy function of the species (in units of J/K).
    std::function<ThermoScalar(double, double)> S;

    /// The apparent standard molar Helmholtz free energy function of the species (in units of J/mol).
    std::function<ThermoScalar(double, double)> A;

    /// The apparent standard molar internal energy function of the species (in units of J/mol).
    std::function<ThermoScalar(double, double)> U;

    /// The apparent standard molar enthalpy function of the species (in units of J/mol).
    std::function<ThermoScalar(double, double)> H;

    /// The apparent standard molar Gibbs free energy function of the species (in units of J/mol).
    std::function<ThermoScalar(double, double)> G;

    /// The standard molar isobaric heat capacity function of the species (in units of J/(mol K)).
    std::function<ThermoScalar(double, double)> Cp;
};

}  // namespace Reaktor


