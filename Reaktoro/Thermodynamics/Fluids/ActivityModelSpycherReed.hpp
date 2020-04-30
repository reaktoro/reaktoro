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

/// The Spycher and Reed (1988) activity model for a gaseous phase.
/// This is an activity model for a gaseous phase supporting only the gases
/// H<sub>2</sub>2O(g), CO<sub>2</sub>(g), and CH<sub>4</sub>(g).
///
/// **References:**
/// - Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2, CH4, H2O
///   and of H2O--CO2--CH4 mixtures: A virial equation treatment for moderate
///   pressures and temperatures applicable to calculations of hydrothermal
///   boiling. Geochimica et Cosmochimica Acta, 52(3), 739�749*.
class ActivityModelSpycherReed : public ActivityModel
{
public:
    /// Construct a default ActivityModelSpycherReed object.
    ActivityModelSpycherReed();

	/// Build the function for activity and thermodynamic excesss property calculations of a phase.
    virtual auto build(const SpeciesList& species) const -> ActivityPropsFn;
};

} // namespace Reaktoro
