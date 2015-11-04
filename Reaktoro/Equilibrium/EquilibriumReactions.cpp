// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "EquilibriumReactions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>

namespace Reaktoro {
namespace {

auto indicesPotentialPrimaryComponents(const Matrix& A) -> Indices
{
    const Index m = A.rows();
    const Index n = A.cols();

    auto nonzeros_in_row = [&](Index irow)
    {
        Index count = 0;
        for(Index j = 0; j < n; ++j)
            if(A(irow, j) != 0) ++count;
        return count;
    };

    auto nonzeros_in_col = [&](Index icol)
    {
        Index count = 0;
        for(Index i = 0; i < m; ++i)
            if(A(i, icol) != 0) ++count;
        return count;
    };

    auto component_weight = [&](Index icol)
    {
        Index weight = 0;
        for(Index i = 0; i < m; ++i)
            if(A(i, icol) != 0) weight += nonzeros_in_row(i);
        return weight;
    };

    auto indices_nonzero_columns_in_row = [&](Index irow)
    {
        Indices indices;
        for(Index j = 0; j < n; ++j)
            if(A(irow, j) != 0) indices.push_back(j);
        return indices;
    };

    auto primary_component_in_row = [&](Index irow)
    {
        Indices nonzero_cols = indices_nonzero_columns_in_row(irow);
        std::sort(nonzero_cols.begin(), nonzero_cols.end(),
            [&](Index l, Index r) {
                if(nonzeros_in_col(l) > nonzeros_in_col(r))
                    return false;
                if(nonzeros_in_col(l) < nonzeros_in_col(r))
                    return true;
                if(component_weight(l) < component_weight(r))
                    return false;
                if(component_weight(l) > component_weight(r))
                    return true;
                if(sum(abs(A.col(l))) < sum(abs(A.col(r))))
                    return true;
                return false;
        });

        return nonzero_cols.front();
    };

    Indices iprimary(m);
    for(Index i = 0; i < m; ++i)
        iprimary[i] = primary_component_in_row(i);

    return iprimary;
}

/// Return a permutation matrix `P` such that the formula matrix
/// `A' = AP` is sorted in order of most potential primary species
auto formulaMatrixPermutation(const Matrix& A) -> PermutationMatrix
{
    const Index m = A.rows();
    const Index n = A.cols();

    Indices indices = range(n);
    Indices iprimary = indicesPotentialPrimaryComponents(A);

    for(Index i = 0; i < m; ++i)
        std::swap(indices[i], iprimary[i]);

    PermutationMatrix P(n);
    std::copy(indices.begin(), indices.end(), P.indices().data());

    return P;
}

} // namespace

struct EquilibriumReactions::Impl
{
    ChemicalSystem system;

    Partition partition;

    DecompositionInfo lu;

    Indices iprimary;

    Indices isecondary;

    std::vector<ReactionEquation> equations;

    Matrix stoichiometric_matrix;
};

EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system)
: EquilibriumReactions(system, Partition(system))
{
}

EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system, const Partition& partition)
{

}

EquilibriumReactions::EquilibriumReactions(const EquilibriumReactions& other)
: pimpl(new Impl(*other.pimpl))
{
}

EquilibriumReactions::~EquilibriumReactions()
{
}

auto EquilibriumReactions::operator=(EquilibriumReactions other) -> EquilibriumReactions&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumReactions::setPrimarySpecies(Indices ispecies) -> void
{

}

auto EquilibriumReactions::setPrimarySpecies(std::vector<std::string> species) -> void
{

}

auto EquilibriumReactions::indicesPrimarySpecies() const -> Indices
{

}

auto EquilibriumReactions::indicesSecondarySpecies() const -> Indices
{

}

auto EquilibriumReactions::equations() const -> std::vector<ReactionEquation>
{

}

auto EquilibriumReactions::stoichiometricMatrix() const -> Matrix
{

}

auto EquilibriumReactions::lu() const -> const DecompositionInfo&
{

}

} // namespace Reaktoro
