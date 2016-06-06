// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "KineticState.hpp"

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct KineticState::Impl
{
};

KineticState::KineticState()
: EquilibriumState(), pimpl(new Impl())
{}

KineticState::KineticState(const ChemicalSystem& system)
: EquilibriumState(system), pimpl(new Impl())
{}

KineticState::KineticState(const ChemicalState& state)
: EquilibriumState(state), pimpl(new Impl())
{}

KineticState::KineticState(const EquilibriumState& state)
: EquilibriumState(state), pimpl(new Impl())
{}

KineticState::KineticState(const KineticState& other)
: EquilibriumState(other), pimpl(new Impl(*other.pimpl))
{}

KineticState::~KineticState()
{}

auto KineticState::operator=(KineticState other) -> KineticState&
{
    ChemicalState::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticState::output(std::string filename) -> void
{
    std::ofstream out(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const KineticState& state) -> std::ostream&
{
    out << static_cast<EquilibriumState>(state);

    return out;
}

} // namespace Reaktoro
