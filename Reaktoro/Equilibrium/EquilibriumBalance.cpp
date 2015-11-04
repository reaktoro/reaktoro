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

#include "EquilibriumBalance.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumReactions.hpp>

namespace Reaktoro {

struct EquilibriumBalance::Impl
{
    /// The system of equilibrium reactions used to regularize the balance matrix
    EquilibriumReactions reactions;

    /// The regularized balance matrix
    Matrix balance_matrix;

    /// The permutation matrix containing the species order corresponding
    /// to the columns of the regularized balance matrix
    PermutationMatrix permutation;

    /// Construct an Impl instance with given reactions
    Impl(const EquilibriumReactions& reactions)
    : reactions(reactions)
    {
        // Initialize the regulaized balance matrix
        const auto& lu = reactions.lu();
        const auto& iprimary = reactions.indicesPrimarySpecies();
        const auto& isecondary = reactions.indicesSecondarySpecies();
        const auto num_primary = iprimary.size();
        const auto num_secondary = isecondary.size();
        const auto num_species = num_primary + num_secondary;
        const auto U1 = lu.U.leftCols(num_primary).triangularView<Eigen::Upper>();
        const auto U2 = lu.U.rightCols(num_secondary);
        balance_matrix.resize(num_primary, num_species);
        balance_matrix.leftCols(num_primary) = identity(num_primary, num_primary);
        balance_matrix.rightCols(num_secondary) = U1.solve(U2);
        cleanRationalNumbers(balance_matrix.data(), balance_matrix.size(), 1e6);

        // Initialize the permutation matrix
        permutation = lu.Q.inverse();
    }
};

EquilibriumBalance::EquilibriumBalance(const ChemicalSystem& system)
: EquilibriumBalance(system, Partition(system))
{
}

EquilibriumBalance::EquilibriumBalance(const ChemicalSystem& system, const Partition& partition)
: EquilibriumBalance(EquilibriumReactions(system, partition))
{
}

EquilibriumBalance::EquilibriumBalance(const EquilibriumReactions& reactions)
: pimpl(new Impl(reactions))
{
}

EquilibriumBalance::EquilibriumBalance(const EquilibriumBalance& other)
: pimpl(new Impl(*other.pimpl))
{
}

EquilibriumBalance::~EquilibriumBalance()
{
}

auto EquilibriumBalance::operator=(EquilibriumBalance other) -> EquilibriumBalance&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumBalance::regularizedMatrix() const -> const Matrix&
{
    return pimpl->balance_matrix;
}

auto EquilibriumBalance::regularizedVector(const Vector& b) const -> Vector
{
    Vector b_bar(b);
    const auto& lu = pimpl->reactions.lu();
    const auto num_primary = pimpl->reactions.indicesPrimarySpecies().size();
    const auto U1 = lu.U.leftCols(num_primary);
    b_bar = lu.P*b_bar;
    b_bar = lu.L.triangularView<Eigen::Lower>().solve(b_bar);
    b_bar = U1.triangularView<Eigen::Upper>().solve(b_bar);
    return b_bar;
}

auto EquilibriumBalance::permutationMatrix() const -> PermutationMatrix
{
    return pimpl->permutation;
}

} // namespace Reaktoro
