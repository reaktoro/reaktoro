// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>

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
    static auto withInitialPriorities(Deque<Index> const& priorities) -> PriorityQueue;

    /// Return a PriorityQueue instance with an initial order and zero priorities.
    static auto withInitialOrder(Deque<Index> const& order) -> PriorityQueue;

    /// Return a PriorityQueue instance with given initial priorities and order.
    /// @note A stable sort algorithm is applied to ensure consistency between
    /// given order and priorities.
    static auto withInitialPrioritiesAndOrder(Deque<Index> const& priorities, Deque<Index> const& order) -> PriorityQueue;

    /// Return the size of the priority queue.
    auto size() const -> Index;

    /// Reset the priorities of the tracked entities.
    auto reset() -> void;

    /// Increment the priority of a tracked entity.
    /// @param identity The index of the tracked entity.
    auto increment(Index identity) -> void;

    /// Extend the queue with the introduction of a new tracked entity.
    auto extend() -> void;

    /// Return the current priorities of each tracked entity in the queue.
    auto priorities() const -> Deque<Index> const&;

    /// Return the current order of the tracked entities in the queue.
    auto order() const -> Deque<Index> const&;

private:
    /// The priorities/usage count of each tracked entity in the priority queue.
    Deque<Index> _priorities;

    /// The order of the tracked entities based on their current priorities.
    Deque<Index> _order;
};

} // namespace Reaktoro
