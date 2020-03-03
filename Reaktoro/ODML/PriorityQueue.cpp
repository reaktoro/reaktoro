// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include "PriorityQueue.hpp"

// C++ includes
#include <algorithm>
#include <numeric>

namespace Reaktoro {

PriorityQueue::PriorityQueue()
{}

auto PriorityQueue::increment(Index ientity) -> void
{
    _rank[ientity] += 1;
    std::sort(_ordering.begin(), _ordering.begin() + ientity + 1,
        [&](Index l, Index r) { return _rank[l] > _rank[r]; });
}

auto PriorityQueue::extend() -> void
{
    _rank.push_back(0);
    _ordering.push_back(_ordering.size());
}

auto PriorityQueue::rank() const -> const std::deque<Index>&
{
    return _rank;
}

auto PriorityQueue::ordering() const -> const std::deque<Index>&
{
    return _ordering;
}

} // namespace Reaktoro
