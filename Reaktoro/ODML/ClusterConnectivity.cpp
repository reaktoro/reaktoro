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

#include "ClusterConnectivity.hpp"

// C++ includes
#include <cassert>
#include <limits>

namespace Reaktoro {

ClusterConnectivity::ClusterConnectivity()
{}

auto ClusterConnectivity::size() const -> Index
{
    return queue.size();
}

auto ClusterConnectivity::extend() -> void
{
    // Extend the priority queue in each row of the connectivity matrix of the clusters
    for(auto& row : matrix)
        row.extend();

    // Extend the priority queue that keeps track the most used clusters
    queue.extend();

    // The ordering of the most used clusters that will be used for the new cluster
    const auto& order = queue.order();

    // Create priorities for the priority queue of the new cluster
    std::deque<Index> priorities(order.size(), 0);

    // Set the priority of the new cluster (the last one) to infinity,
    // to ensure it is always the first one to be visited as we
    // go from cluster to cluster
    priorities.back() = std::numeric_limits<Index>::infinity();

    // Append the new priority queue for the new cluster
    matrix.push_back(PriorityQueue::withInitialPrioritiesAndOrder(priorities, order));
}

auto ClusterConnectivity::increment(Index icluster, Index jcluster) -> void
{
    // Only jcluster needs to be bounded, because icluster >= size() has a specific logic
    assert(jcluster < size());

    // Increment jcluster when starting from icluster (if icluster is below number of clusters!)
    if(icluster < size())
        matrix[icluster].increment(jcluster);

    // Increment usage count of jcluster
    queue.increment(jcluster);
}

auto ClusterConnectivity::order(Index icluster) const -> const std::deque<Index>&
{
    return icluster < size() ? matrix[icluster].order() : queue.order();
}

} // namespace Reaktoro

