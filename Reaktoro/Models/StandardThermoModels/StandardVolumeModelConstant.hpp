// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/StandardVolumeModel.hpp>

namespace Reaktoro {

/// The parameters in the constant model for calculating standard volume of a product species in a formation reaction.
struct StandardVolumeModelParamsConstant
{
    /// The constant standard molar volume @f$V^{\circ}@f$ of the product species (in m3/mol).
    Param V0;
};

/// Return a function that calculates the standard volume of a product species using a constant model.
auto StandardVolumeModelConstant(const StandardVolumeModelParamsConstant& params) -> StandardVolumeModel;

} // namespace Reaktoro
