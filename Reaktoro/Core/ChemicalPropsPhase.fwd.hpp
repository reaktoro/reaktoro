// Reaktoro is a unified framework for modeling chemically reactive phases.
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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The base type for primary chemical property data of a phase from which others are computed.
template<typename Real, typename Array>
struct ChemicalPropsPhaseBaseData;

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseData = ChemicalPropsPhaseBaseData<real, ArrayXr>;

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseDataRef = ChemicalPropsPhaseBaseData<real&, ArrayXrRef>;

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseDataConstRef = ChemicalPropsPhaseBaseData<const real&, ArrayXrConstRef>;


/// The base type for chemical properties of a phase and its species.
template<typename R, typename A>
class ChemicalPropsPhaseBase;

/// The chemical properties of a phase and its species.
using ChemicalPropsPhase = ChemicalPropsPhaseBase<real, ArrayXr>;

/// The non-const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseRef = ChemicalPropsPhaseBase<real&, ArrayXrRef>;

/// The const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseConstRef = ChemicalPropsPhaseBase<const real&, ArrayXrConstRef>;

} // namespace Reaktoro
