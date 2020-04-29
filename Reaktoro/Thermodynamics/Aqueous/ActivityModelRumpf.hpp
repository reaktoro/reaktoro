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

/// The activity model of Duan and Sun (2003) for a dissolved gas.
/// **References:**
///  - Rumpf, B., Nicolaisen, H., Ocal, C., & Maurer, G. (1994). Solubility of
///    carbon dioxide in aqueous mixtures of sodium chloride: Experimental
///    results and correlation. Journal of Solution Chemistry, 23(3), 431–448*.
class ActivityModelRumpf : public ActivityModel
{
public:
    /// Construct a default ActivityModelRumpf object for dissolved gas CO<sub>2</sub>(aq).
    ActivityModelRumpf();

    /// Construct a ActivityModelRumpf object with given dissolved gas formula.
    ActivityModelRumpf(String gas);

    /// Build the function for activity and thermodynamic excesss property calculations of a phase.
    virtual auto build(const SpeciesList& species) const -> ActivityPropsFn;

private:
    /// The chemical formula of the dissolved gas.
    String gas = "CO2";
};

} // namespace Reaktoro
