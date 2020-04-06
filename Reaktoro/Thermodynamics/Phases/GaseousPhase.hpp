// Reaktoro is a unified framework for modeling chemically reactive systems.
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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Phases/FluidPhase.hpp>

namespace Reaktoro {

using GeneralMixture = GeneralMixture;

class GaseousPhase
    : public FluidPhase
{
public:
    /// Construct a default GaseousPhase instance.
    GaseousPhase()
        : FluidPhase("Gaseous", StateOfMatter::Gas)
    {}

    /// Construct a GaseousPhase instance with given gaseous mixture.
    /// The Peng-Robinson equation of state is chosen by default to calculate the
    /// thermodynamic and chemical properties of this GaseousPhase object.
    explicit GaseousPhase(const GeneralMixture& mixture)
        : FluidPhase(mixture, "Gaseous", StateOfMatter::Gas)
    {}
};

} // namespace Reaktoro
