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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

/// The parameters in the Extended UNIQUAC activity model for aqueous electrolyte solutions.
struct ActivityModelParamsExtendedUNIQUAC
{
    /// The volume parameters \eq{r_i} for the aqueous species.
    Tuples<String, Vec<real>> r;

    /// The surface parameters \eq{r_i} for the aqueous species.
    Tuples<String, Vec<real>> q;

    /// The binary interaction parameters \eq{u_{ij}} among the aqueous species.
    Tuples<String, String, Vec<real>> u;
};

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC() -> ActivityModelGenerator;

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC(ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModelGenerator;

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC(Params const& params) -> ActivityModelGenerator;

} // namespace Reaktoro
