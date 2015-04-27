// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "ChemicalPlot.hpp"

// C++ includes
#include <iomanip>

// Reaktoro includes
#include <Reaktoro/Common/Gnuplot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

struct ChemicalPlot::Impl
{
    ChemicalSystem system;

    ReactionSystem reactions;

    Gnuplot gnuplot;

    std::stringstream data;
//
//    std::ofstream file;

    Impl()
    {}

    Impl(const ChemicalSystem& system)
    : system(system)
    {}

    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions)
    {}
};

ChemicalPlot::ChemicalPlot()
: pimpl(new Impl())
{}

ChemicalPlot::ChemicalPlot(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalPlot::ChemicalPlot(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalPlot::~ChemicalPlot()
{}

auto ChemicalPlot::operator<<(std::string str) -> ChemicalPlot&
{
    pimpl->gnuplot << str;
    return *this;
}

auto ChemicalPlot::operator<<(const std::stringstream& ss) -> ChemicalPlot&
{
    pimpl->gnuplot << ss;
    return *this;
}

auto ChemicalPlot::update(double t, const ChemicalState& state) -> void
{
    auto& data = pimpl->data;
    data << std::left << std::setw(20) << t;
    for(auto quantity : what)
        data << std::left << std::setw(20) << extract(state, quantity);
    data << std::endl;
}

} // namespace Reaktoro
