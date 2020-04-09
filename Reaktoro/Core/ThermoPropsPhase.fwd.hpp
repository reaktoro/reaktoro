// Reaktoro is a unified framework for modeling chemically reactive phases.
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

/// The base type for primary standard thermodynamic property data of a phase.
template<typename Real, typename Array>
struct ThermoPropsPhaseBaseData;

/// The primary standard thermodynamic property data of a phase.
using ThermoPropsPhaseData = ThermoPropsPhaseBaseData<real, ArrayXr>;

/// The primary standard thermodynamic property data of a phase.
using ThermoPropsPhaseDataRef = ThermoPropsPhaseBaseData<real&, ArrayXrRef>;

/// The primary standard thermodynamic property data of a phase.
using ThermoPropsPhaseDataConstRef = ThermoPropsPhaseBaseData<const real&, ArrayXrConstRef>;


/// The base type for standard thermodynamic properties of a phase and its species.
template<typename R, typename A>
class ThermoPropsPhaseBase;

/// The standard thermodynamic properties of a phase and its species.
using ThermoPropsPhase = ThermoPropsPhaseBase<real, ArrayXr>;

/// The non-const view to the standard thermodynamic properties of a phase and its species.
using ThermoPropsPhaseRef = ThermoPropsPhaseBase<real&, ArrayXrRef>;

/// The const view to the standard thermodynamic properties of a phase and its species.
using ThermoPropsPhaseConstRef = ThermoPropsPhaseBase<const real&, ArrayXrConstRef>;

} // namespace Reaktoro
