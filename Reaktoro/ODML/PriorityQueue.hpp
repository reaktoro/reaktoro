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

// A queue organized based on priorities that can change dynamically.
class PriorityQueue
{
public:
    /// Construct a default instance of PriorityQueue.
    PriorityQueue();

    /// Return a PriorityQueue instance with given initial size.
    static auto withInitialSize(Index size) -> PriorityQueue;

    /// Return a PriorityQueue instance with given initial priorities.
    static auto withInitialPriorities(const std::deque<Index>& priorities) -> PriorityQueue;

    /// Return a PriorityQueue instance with an initial order and zero priorities.
    static auto withInitialOrder(const std::deque<Index>& order) -> PriorityQueue;

    /// Return a PriorityQueue instance with given initial priorities and order.
    /// @note A stable sort algorithm is applied to ensure consistency between
    /// given order and priorities.
    static auto withInitialPrioritiesAndOrder(const std::deque<Index>& priorities, const std::deque<Index>& order) -> PriorityQueue;

    /// Return the size of the priority queue.
    auto size() const -> Index;

    /// Reset the priorities of the tracked entities.
    auto reset() -> void;

    /// Increment the priority of a tracked entity.
    /// @param ientity The index of the tracked entity.
    auto increment(Index ientity) -> void;

    /// Extend the queue with the introduction of a new tracked entity.
    auto extend() -> void;

    /// Return the current priorities of each tracked entity in the queue.
    auto priorities() const -> const std::deque<Index>&;

    /// Return the current order of the tracked entities in the queue.
    auto order() const -> const std::deque<Index>&;

private:
    /// The priorities/usage count of each tracked entity in the priority queue.
    std::deque<Index> _priorities;

    /// The order of the tracked entities based on their current priorities.
    std::deque<Index> _order;
};

} // namespace Reaktoro
