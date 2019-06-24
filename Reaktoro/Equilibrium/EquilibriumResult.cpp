// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "EquilibriumResult.hpp"

namespace Reaktoro {

auto SmartEquilibriumResult::EstimateStatistics::operator+=(const SmartEquilibriumResult::EstimateStatistics& other)
-> SmartEquilibriumResult::EstimateStatistics&
{
    time_estimate += other.time_estimate;
    time_search += other.time_search;
    time_mat_vect_mult += other.time_mat_vect_mult;
    time_acceptance += other.time_acceptance;
    return *this;
}
auto SmartEquilibriumResult::LearnStatistics::operator+=(const SmartEquilibriumResult::LearnStatistics& other)
        -> SmartEquilibriumResult::LearnStatistics&
{
    time_learn += other.time_learn;
    time_store += other.time_store;
    time_gibbs_min += other.time_gibbs_min;
    return *this;
}

auto SmartEquilibriumResult::operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&
{
    succeeded       = other.succeeded;
    learn_stats     += other.learn_stats;
    estimate_stats  += other.estimate_stats;
    return *this;
}

auto SmartEquilibriumResult::addLearningIndex(const Index & index) -> void{
    learning_states_indx.emplace_back(index);
}

auto EquilibriumResult::operator+=(const EquilibriumResult& other) -> EquilibriumResult&
{
    optimum += other.optimum;
    smart   += other.smart;
    stats   += other.stats;
    return *this;
}

} // namespace Reaktoro
