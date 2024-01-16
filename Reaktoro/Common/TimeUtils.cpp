// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "TimeUtils.hpp"

namespace Reaktoro {

auto time() -> Time
{
    return std::chrono::high_resolution_clock::now();
}

auto elapsed(Time const& end, Time const& begin) -> double
{
    return std::chrono::duration<double>(end - begin).count();
}

auto elapsed(Time const& begin) -> double
{
    return elapsed(time(), begin);
}

Stopwatch::Stopwatch()
{}

auto Stopwatch::start() -> void
{
    mstart = Reaktoro::time();
}

auto Stopwatch::pause() -> void
{
    melapsed += elapsed(mstart);
}

auto Stopwatch::reset() -> void
{
    mstart = {};
    melapsed = {};
}

auto Stopwatch::time() const -> double
{
    return melapsed;
}

} // namespace Reaktoro
