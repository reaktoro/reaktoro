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
    time_search += other.time_search;
    time_matrix_vector_mult += other.time_matrix_vector_mult;
    time_acceptance_test += other.time_acceptance_test;
    return *this;
}
auto SmartEquilibriumResult::LearnStatistics::operator+=(const SmartEquilibriumResult::LearnStatistics& other)
        -> SmartEquilibriumResult::LearnStatistics&
{
    time_store += other.time_store;
    time_gibbs_minimization += other.time_gibbs_minimization;
    return *this;
}

auto SmartEquilibriumResult::operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&
{
    learn_stats     += other.learn_stats;
    estimate_stats  += other.estimate_stats;
    return *this;
}
auto EquilibriumResult::operator+=(const EquilibriumResult& other) -> EquilibriumResult&
{
    optimum += other.optimum;
    smart   += other.smart;
    return *this;
}

} // namespace Reaktoro
