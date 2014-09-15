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

namespace Reaktor {

/// A type used to store 
/// @see ThermoState
struct ThermoStateInterpolator
{
    /// The interpolator of the standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
    BilinearInterpolator volume;

    /// The interpolator of the standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
    BilinearInterpolator entropy;

    /// The interpolator of the apparent standard molar Helmholtz free energy of formation @f$\Delta A_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator helmholtz;

    /// The interpolator of the apparent standard molar internal energy of formation @f$\Delta U_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator internal_energy;

    /// The interpolator of the apparent standard molar enthalpy of formation @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator enthalpy;

    /// The interpolator of the apparent standard molar Gibbs free energy of formation @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
    BilinearInterpolator gibbs;

    /// The interpolator of the standard molar isobaric heat capacity @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol K))
    BilinearInterpolator cp;
};

}  // namespace Reaktor


