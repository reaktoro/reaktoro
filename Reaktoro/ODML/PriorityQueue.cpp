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

#include "PriorityQueue.hpp"

// C++ includes
#include <algorithm>
#include <cassert>
#include <numeric>

namespace Reaktoro {

PriorityQueue::PriorityQueue()
{}

auto PriorityQueue::withInitialSize(Index size) -> PriorityQueue
{
    PriorityQueue queue;
    queue._priorities.resize(size, 0);
    queue._order.resize(size);
    std::iota(queue._order.begin(), queue._order.end(), 0);
    return queue;
}

auto PriorityQueue::withInitialPriorities(Deque<Index> const& priorities) -> PriorityQueue
{
    const auto size = priorities.size();
    PriorityQueue queue = PriorityQueue::withInitialSize(size);
    queue._priorities = priorities;
    std::sort(queue._order.begin(), queue._order.end(),
        [&](Index l, Index r) { return priorities[l] > priorities[r]; });
    return queue;
}

auto PriorityQueue::withInitialOrder(Deque<Index> const& order) -> PriorityQueue
{
    const auto size = order.size();
    PriorityQueue queue;
    queue._priorities.resize(size, 0);
    queue._order = order;
    return queue;
}

auto PriorityQueue::withInitialPrioritiesAndOrder(Deque<Index> const& priorities, Deque<Index> const& order) -> PriorityQueue
{
    assert(priorities.size() == order.size());
    PriorityQueue queue;
    queue._priorities = priorities;
    queue._order = order;
    std::stable_sort(queue._order.begin(), queue._order.end(),
        [&](Index l, Index r) { return priorities[l] > priorities[r]; });
    return queue;
}

auto PriorityQueue::size() const -> Index
{
    return _priorities.size();
}

auto PriorityQueue::reset() -> void
{
    std::fill(_priorities.begin(), _priorities.end(), 0);
    std::iota(_order.begin(), _order.end(), 0);
}

auto PriorityQueue::increment(Index identity) -> void
{
    // == EXAMPLE OF WHAT HAPPENS IN THIS METHOD ==
    // PRIORITIES BEFORE INCREMENTING: 13  5  3  2  2 (2) 1  --- incrementing from 2 to 3
    //  PRIORITIES AFTER INCREMENTING: 13  5  3 [2] 2 (3) 1  --- (3) needs to be swapped with [2]
    //       PRIORITIES AFTER SORTING: 13  5  3 (3) 2 [2] 1
    _priorities[identity] += 1;

    // == CORNER CASES:
    // IF identity = 0 (the priority of the very first element is increased)
    // NO SORTING IS APPLIED (THE CHANGE IN PRIORITIES DIDN'T AFFECT THE ORDER)
    // PRIORITIES BEFORE INCREMENTING: (13)  5  3  2  2  2  1    --- incrementing from 13 to 14
    //  PRIORITIES AFTER INCREMENTING: (14)  5  3  2  2  3  1  --- (14) does not need to be swapped with anyone

    // == CORNER CASES:
    // IF _priorities[identity - 1] >= _priorities[identity]
    // NO SORTING IS APPLIED (THE CHANGE IN PRIORITIES DIDN'T AFFECT THE ORDER)
    // PRIORITIES BEFORE INCREMENTING: 13  5 (3)  2  2  2  1    --- incrementing from 3 to 4
    //  PRIORITIES AFTER INCREMENTING: 13  5 (4)  2  2  3  1  --- (4) is still smaller than 5, it does not need to be swapped with anyone

    // Avoid sorting if the identity == 0 (the priority of the very first element is increased)
    // or if the increase of priority does not ruin the order in _priorities
    if(!((identity == 0) || (_priorities[identity - 1] >= _priorities[identity]))){
        std::sort(_order.begin(), _order.begin() + identity + 1,
                  [&](Index l, Index r) { return _priorities[l] > _priorities[r]; });
    }
    //std::sort(_order.begin(), _order.begin() + identity + 1,
    //    [&](Index l, Index r) { return _priorities[l] > _priorities[r]; });
}

auto PriorityQueue::extend() -> void
{
    _priorities.push_back(0);
    _order.push_back(_order.size());
}

auto PriorityQueue::priorities() const -> Deque<Index> const&
{
    return _priorities;
}

auto PriorityQueue::order() const -> Deque<Index> const&
{
    return _order;
}

} // namespace Reaktoro
