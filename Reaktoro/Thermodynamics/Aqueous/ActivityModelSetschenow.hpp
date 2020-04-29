// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

namespace Reaktoro {

/// The Setschenow activity model for a neutral aqueous species.
class ActivityModelSetschenow : public ActivityModel
{
public:
    /// Construct a ActivityModelSetschenow object.
    /// @param neutral The formula of the neutral aqueous species (e.g., `NaCl`).
    /// @param b The Setschenow *b* coefficient.
    ActivityModelSetschenow(String neutral, real b);

	/// Build the function for activity and thermodynamic excesss property calculations of a phase.
    virtual auto build(const SpeciesList& species) const -> ActivityPropsFn;

private:
    /// The formula of the neutral aqueous species.
    String neutral;

    /// The Setschenow *b* coefficient.
    real b;
};


} // namespace Reaktoro
