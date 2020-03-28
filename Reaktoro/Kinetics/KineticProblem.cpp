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

#include "KineticProblem.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

struct KineticProblem::Impl
{
    ReactionSystem reactions;

    ChemicalSystem system;

    Partition partition;

    double temperature;

    double pressure;

    VectorXr n;

    double initial_time;

    double final_time;

    Impl()
    {}

    Impl(const ReactionSystem& reactions)
    : reactions(reactions), system(reactions.system()),
      temperature(INFINITY), pressure(INFINITY),
      initial_time(0), final_time(INFINITY)
    {}
};

KineticProblem::KineticProblem(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

KineticProblem::KineticProblem(const KineticProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticProblem::~KineticProblem()
{}

auto KineticProblem::operator=(KineticProblem other) -> KineticProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticProblem::setTemperature(double val) -> KineticProblem&
{
    pimpl->temperature = val;
    return *this;
}

auto KineticProblem::setPressure(double val) -> KineticProblem&
{
    pimpl->pressure = val;
    return *this;
}

auto KineticProblem::temperature() const -> double
{
    return pimpl->temperature;
}

auto KineticProblem::pressure() const -> double
{
    return pimpl->pressure;
}

} // namespace Reaktoro

