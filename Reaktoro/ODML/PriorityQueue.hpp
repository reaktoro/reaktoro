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

#pragma once

// C++ includes
#include <deque>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// A priority queue data structure based on usage count.
class PriorityQueue
{
public:
    /// Construct a default instance of PriorityQueue.
    PriorityQueue();

    /// Increment the rank of a tracked entity to increase its priority.
    /// @param ientity The index of the tracked entity.
    auto increment(Index ientity) -> void;

    /// Extend the priority queue with the introduction of a new tracked entity.
    auto extend() -> void;

    /// Return the current rank of each tracked entity in the priority queue.
    auto rank() const -> const std::deque<Index>&;

    /// Return the current ordering of the tracked entities in the priority queue based on their rank.
    auto ordering() const -> const std::deque<Index>&;

private:
    /// The rank/usage count of each tracked entity in the priority queue.
    std::deque<Index> _rank;

    /// The ordering of the tracked entities based on their current rank.
    std::deque<Index> _ordering;
};

} // namespace Reaktoro
