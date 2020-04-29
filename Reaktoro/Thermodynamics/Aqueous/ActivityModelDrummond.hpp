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

/// The activity model of Drummond (1981) for dissolved gases.
/// In the activity model of Drummond (1981), the activity coefficient of a
/// dissolved gas (e.g., CO<sub>2</sub>(aq), O<sub>2</sub>(aq),
/// H<sub>2</sub>S(aq)) is computed using:
/// @eqc{\ln\gamma_{i}=\left(a_{1}+a_{2}T+\frac{a_{3}}{T}\right)-(a_{4}+a_{5}T)\frac{I}{I+1},}
/// where I is the ionic strength of the aqueous solution (in molal) and T is
/// temperature (in K).
///
/// **Reference:**
///
/// *Drummond, S. E. (1981). Boiling and mixing of hydrothermal fluids:
/// chemical effects on mineral precipitation. Pennsylvania State University.*
class ActivityModelDrummond
{
public:
    /// The parameters in the Drummond (1981) activity model for dissolved gases.
    /// The default values correspond to dissolved gas CO<sub>2</sub>(aq).
    struct Params
    {
        real a1 = -1.0312;    ///< The coefficient @eq{a_1} in the activity model.
        real a2 =  1.2806e-3; ///< The coefficient @eq{a_2} in the activity model.
        real a3 =  255.9;     ///< The coefficient @eq{a_3} in the activity model.
        real a4 =  0.4445;    ///< The coefficient @eq{a_4} in the activity model.
        real a5 = -1.6060e-3; ///< The coefficient @eq{a_5} in the activity model.
    };

    /// Construct a default ActivityModelDrummond object for dissolved gas CO<sub>2</sub>(aq).
    ActivityModelDrummond();

    /// Construct a ActivityModelDrummond object with given dissolved gas formula.
    ActivityModelDrummond(String gas);

    /// Construct a ActivityModelDrummond object with given dissolved gas formula and model parameters.
    ActivityModelDrummond(String gas, Params params);

    /// Construct the activity function with given list of species.
    auto operator()(const SpeciesList& species) const -> ActivityPropsFn;

private:
    /// The chemical formula of the dissolved gas.
    String gas = "CO2";

    /// The parameters of the activity model.
    Params params;
};

} // namespace Reaktoro
