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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return an activity model that applies an ionic strengh and pressure
/// correction on the activity coefficients of aqueous solutes to produce
/// consistent results with PHREEQC. This is needed because PHREEQC introduces a
/// dependency on ionic strength for the equilibrium constants of reactions and
/// this is incompatible with the design decision in Reaktoro to make these
/// quantities a function of temperature and presssure, not composition. Thus,
/// the ionic strenth correction used in PHREEQC needs to be transfer in
/// Reaktoro to the activity coefficients of the aqueous solutes. Thus, the use
/// of this activity model needs to be done via a chain of activity models. This
/// corrective activity model should be the last one in the chainned model. Example:
/// ~~~c++
/// AqueousPhase solution(speciate("Na Cl C"));
/// solution.setActivityModel(chain(
///     ActivityModelPitzer(),
///     ActivityModelPhreeqcIonicStrengthPressureCorrection(),
/// ));
/// ~~~
auto ActivityModelPhreeqcIonicStrengthPressureCorrection() -> ActivityModelGenerator;

} // namespace Reaktoro
