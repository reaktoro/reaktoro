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
#include <algorithm>
#include <cassert>
#include <numeric>

namespace Reaktoro {

ClusterConnectivity::ClusterConnectivity()
{}

auto ClusterConnectivity::size() const -> Index
{
    return _rank.size();
}

auto ClusterConnectivity::extend() -> void
{
    // The current number of rows and columns in the matrix
    const auto size = _rank.size();

    // Append a zero count at the end of each rank row for the new cluster (not yet used)
    for(auto& row : _rank)
        row.push_back(0);

    // Append the index of new cluster to the end of each ordering row (it should be the last one tried, as it is new)
    for(auto& row : _ordering)
        row.push_back(size);

    // Create a new rank row (with zeros)
    std::deque<Index> new_rank_row(size + 1, 0);

    // Set the last entry in the last counting row to infinity, to give it first priority always
    new_rank_row.back() = std::numeric_limits<Index>::infinity();

    // Create a new ordering row
    std::deque<Index> new_ordering_row(size + 1);

    // Set first index to the index of current cluster created
    new_ordering_row[0] = size;

    // Set the indices of the other clusters in ascending order, starting from 0
    std::iota(new_ordering_row.begin() + 1, new_ordering_row.end(), 0);

    // As initial ordering, consider the total usage count of clusters beyond the first
    std::sort(new_ordering_row.begin() + 1, new_ordering_row.end(),
        [&](Index l, Index r) { return _rank_total[l] > _rank_total[r]; });

    // Push back the new row of rank and ordering
    _rank.push_back(new_rank_row);
    _ordering.push_back(new_ordering_row);

    // Update the total rank/usage count list and the ordering based on total usage counts
    _rank_total.push_back(0);
    _ordering_total.push_back(size);
}

auto ClusterConnectivity::increment(Index icluster, Index jcluster) -> void
{
    // Only jcluster needs to be bounded, because icluster >= size() has a specific logic
    assert(jcluster < size());

    // Check if index of starting cluster is below number of clusters
    if(icluster < size())
    {
        // Increment rank/usage count of connectivity from icluster to jcluster
        _rank[icluster][jcluster] += 1;

        // Update ordering of clusters for the i-th cluster to reflect change in rank of j-th cluster
        std::sort(_ordering[icluster].begin(), _ordering[icluster].begin() + jcluster + 1); // TODO This sort operation can be made faster, since the deque was already sorted
    }

    // Increment total rank/usage count of cluster jcluster
    _rank_total[jcluster] += 1;

    // Update ordering of clusters based on their total usage counts
    std::sort(_ordering_total.begin(), _ordering_total.begin() + jcluster + 1); // TODO This sort operation can be made faster, since the deque was already sorted
}

auto ClusterConnectivity::ordering(Index icluster) const -> const std::deque<Index>&
{
    return icluster < size() ? _ordering[icluster] : _ordering_total;
}

} // namespace Reaktoro

