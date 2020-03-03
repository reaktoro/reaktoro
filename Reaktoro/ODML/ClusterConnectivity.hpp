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

// The connectivity matrix of the clusters.
class ClusterConnectivity
{
public:
    /// Construct a default instance of ClusterConnectivity.
    ClusterConnectivity();

    /// Return number of currently tracked clusters.
    auto size() const -> Index;

    /// Extend the connectivity matrix following creation of a new cluster.
    auto extend() -> void;

    /// Increment the rank/usage count for the connectivity from one cluster to another.
    /// -----------------------------------------------------------------------
    /// @param icluster The index of the starting cluster.
    /// @param jcluster The index of the cluster which usage count is incremented.
    /// -----------------------------------------------------------------------
    /// @note If index `icluster` is equal or greater than number of clusters,
    /// then usage count for `jcluster` is still incremented, but in a separate
    /// list, which can be used when a starting cluster cannot be specified and
    /// ordering of clusters solely based on their usage counts is used instead.
    auto increment(Index icluster, Index jcluster) -> void;

    /// Return the ordering of clusters for a given starting cluster.
    /// -----------------------------------------------------------------------
    /// @param icluster The index of the starting cluster.
    /// -----------------------------------------------------------------------
    /// @note If index `icluster` is equal or greater than number of clusters,
    /// then an ordering based on usage count of clusters is returned.
    auto ordering(Index icluster) const -> const std::deque<Index>&;

private:
    /// The rank/usage count of each connectivity pair of clusters.
    std::deque<std::deque<Index>> _rank;

    /// The ordering of clusters for each cluster based on their connectivity usage counts.
    std::deque<std::deque<Index>> _ordering;

    /// The total rank/usage count of each cluster.
    std::deque<Index> _rank_total;

    /// The ordering of clusters for each cluster based on their total usage counts.
    std::deque<Index> _ordering_total;
};

} // namespace Reaktoro
