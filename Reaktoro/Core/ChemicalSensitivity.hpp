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
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// The type that stores sensitivity matrices of a chemical state.
class ChemicalSensitivity
{
public:
    /// The sensitivity of the molar amounts of species w.r.t. temperature (in units of mol/K).
    Vector dndT;

    /// The sensitivity of the molar amounts of species w.r.t. pressure (in units of mol/Pa).
    Vector dndP;

    /// The sensitivity of the molar amounts of species w.r.t. molar amounts of elements (in units of mol/mol).
    Matrix dndb;

    /// The sensitivity of the molar amounts of species w.r.t. time (in units of mol/s).
    Vector dndt;
};

} // namespace Reaktoro
