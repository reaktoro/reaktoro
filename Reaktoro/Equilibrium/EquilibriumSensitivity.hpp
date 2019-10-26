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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A type that contains the sensitivity data of the equilibrium state.
/// The sensitivity of the equilibrium state contains derivatives of the species
/// amounts with respect to temperature, pressure, and amounts of the elements.
/// It is an important information for implicit numerical methods, since these
/// derivatives allow the calculation to converge faster to the solution.
struct EquilibriumSensitivity
{
    /// The partial derivatives @f$\left.\frac{\partial n}{\partial T}\right|_{P,b}@f$ (in units of mol/K).
    Vector dndT;

    /// The partial derivatives @f$\left.\frac{\partial y/RT}{\partial T}\right|_{P,b}@f$ (in units of 1/K).
    Vector dydT;

    /// The partial derivatives @f$\left.\frac{\partial z/RT}{\partial T}\right|_{P,b}@f$ (in units of 1/K).
    Vector dzdT;

    /// The partial derivatives @f$\left.\frac{\partial r/RT}{\partial T}\right|_{P,b}@f$ (in units of 1/K).
    Vector drdT;

    /// The partial derivatives @f$\left.\frac{\partial n}{\partial P}\right|_{T,b}@f$ (in units of mol/Pa).
    Vector dndP;

    /// The partial derivatives @f$\left.\frac{\partial y/RT}{\partial P}\right|_{T,b}@f$ (in units of 1/Pa).
    Vector dydP;

    /// The partial derivatives @f$\left.\frac{\partial z/RT}{\partial P}\right|_{T,b}@f$ (in units of 1/Pa).
    Vector dzdP;

    /// The partial derivatives @f$\left.\frac{\partial r/RT}{\partial P}\right|_{T,b}@f$ (in units of 1/Pa).
    Vector drdP;

    /// The partial derivatives @f$\left.\frac{\partial n}{\partial b}\right|_{T, P}@f$ (in units of mol/mol).
    Matrix dndb;

    /// The partial derivatives @f$\left.\frac{\partial y/RT}{\partial b}\right|_{T, P}@f$ (in units of 1/mol).
    Matrix dydb;

    /// The partial derivatives @f$\left.\frac{\partial z/RT}{\partial b}\right|_{T, P}@f$ (in units of 1/mol).
    Matrix dzdb;

    /// The partial derivatives @f$\left.\frac{\partial r/RT}{\partial b}\right|_{T, P}@f$ (in units of 1/mol).
    Matrix drdb;
};

} // namespace Reaktoro
