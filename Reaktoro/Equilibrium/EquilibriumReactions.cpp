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
#include <Reaktoro/Common/Exception.hpp>
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

/// Return a permutation matrix `P` such that the columns of the formula
/// matrix `A' = AP` is reordered where the first columns are swaped with
/// the ones in `ipriority`.
auto formulaMatrixPermutation(const Matrix& A, const Indices& ipriority) -> PermutationMatrix
{
    const Index m = A.rows();
    const Index n = A.cols();

    Indices indices = range(n);

    for(Index i = 0; i < m; ++i)
        std::swap(indices[i], indices[ipriority[i]]);

    PermutationMatrix permutation(n);
    std::copy(indices.begin(), indices.end(), permutation.indices().data());

    return permutation;
}

/// Decompose the formula matrix `A` in LU factors `PAQ = LU`.
/// The columns of the formula matrix `A` is first reordered where
/// the column indices given in `ipriority` are swaped with the fist
/// columns of `A`.
auto luFormulaMatrix(Matrix A, const Indices& ipriority) -> DecompositionLU
{
    PermutationMatrix W = formulaMatrixPermutation(A, ipriority);
    A = A*W; // reorder the columns of A with column indices in ipriority coming first
    DecompositionLU lu = Reaktoro::lu(A); // perform the LU decomposition using partial (row-wise only) pivoting
    lu.Q = W*lu.Q; // combine the permutation matrix Q from the LU decomposition with W
    return lu;
}

} // namespace

struct EquilibriumReactions::Impl
{
    ChemicalSystem system;

    Partition partition;

    Matrix Ae;

    Indices iequilibrium;

    DecompositionLU lu;

    Indices iprimary;

    Indices isecondary;

    Matrix stoichiometric_matrix;

    std::vector<ReactionEquation> equations;

    Impl(const ChemicalSystem& system)
    : Impl(system, Partition(system))
    {}

    Impl(const ChemicalSystem& system, const Partition& partition)
    : system(system), partition(partition)
    {
        Ae = partition.formulaMatrixEquilibriumSpecies();
        iequilibrium = partition.indicesEquilibriumSpecies();
        Indices ipriority = indicesPotentialPrimaryComponents(Ae);
        initialize(ipriority);
    }

    auto initialize(const Indices& ipriority) -> void
    {
        // The number of species in the equilibrium partition
        const Index num_species = Ae.cols();

        // Perform the LU decomposition of the formula matrix of the equilibrium partition
        lu = luFormulaMatrix(Ae, ipriority);

        // The number of primary and secondary species
        const Index num_primary = lu.rank;
        const Index num_secondary = num_species - lu.rank;

        // Initialize global indices of primary and secondary species in the equilibrium partition
        iprimary.resize(num_primary);
        isecondary.resize(num_secondary);
        for(Index i = 0; i < num_primary; ++i)
            iprimary[i] = iequilibrium[lu.Q.indices()[i]];
        for(Index i = 0; i < num_secondary; ++i)
            isecondary[i] = iequilibrium[lu.Q.indices()[i + num_primary]];

        // Initialize the stoichiometric matrix
        const auto U1 = lu.U.leftCols(num_primary).triangularView<Eigen::Upper>();
        const auto U2 = lu.U.rightCols(num_secondary);
        stoichiometric_matrix.resize(num_secondary, num_species);
        stoichiometric_matrix.leftCols(num_primary) = tr(U1.solve(U2));
        stoichiometric_matrix.rightCols(num_secondary) = -identity(num_secondary, num_secondary);
        stoichiometric_matrix = stoichiometric_matrix * lu.Q.inverse();

        // The maximum allowed denominator in the rational numbers forming the stoichiometric matrix
        const long maxden = 1e6;

        // Clean the stoichiometric matrix and the LU factors from round-off errors
        cleanRationalNumbers(stoichiometric_matrix.data(), stoichiometric_matrix.size(), maxden);
        cleanRationalNumbers(lu.L.data(), lu.L.size(), maxden);
        cleanRationalNumbers(lu.U.data(), lu.U.size(), maxden);

        // Initialize the system of cannonical equilibrium reactions
        equations.clear(); equations.reserve(num_secondary);
        for(Index i = 0; i < num_secondary; ++i)
        {
            std::map<std::string, double> equation;
            for(Index j = 0; j < num_species; ++j)
                if(stoichiometric_matrix(i, j) != 0)
                    equation.emplace(system.species(iequilibrium[j]).name(), stoichiometric_matrix(i, j));
            equations.push_back(ReactionEquation(equation));
        }
    }
};

EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system)
: EquilibriumReactions(system, Partition(system))
{
}

EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system, const Partition& partition)
: pimpl(new Impl(system, partition))
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

auto EquilibriumReactions::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumReactions::partition() const -> const Partition&
{
    return pimpl->partition;
}


auto EquilibriumReactions::setPrimarySpecies(Indices ispecies) -> void
{
    pimpl->initialize(ispecies);
}

auto EquilibriumReactions::setPrimarySpecies(std::vector<std::string> species) -> void
{
    Indices ispecies = pimpl->system.indicesSpecies(species);
    setPrimarySpecies(ispecies);
}

auto EquilibriumReactions::indicesPrimarySpecies() const -> Indices
{
    return pimpl->iprimary;
}

auto EquilibriumReactions::indicesSecondarySpecies() const -> Indices
{
    return pimpl->isecondary;
}

auto EquilibriumReactions::equations() const -> std::vector<ReactionEquation>
{
    return pimpl->equations;
}

auto EquilibriumReactions::stoichiometricMatrix() const -> Matrix
{
    return pimpl->stoichiometric_matrix;
}

auto EquilibriumReactions::lu() const -> const DecompositionLU&
{
    return pimpl->lu;
}

} // namespace Reaktoro
