// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "EquilibriumInverseProblem.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

struct EquilibriumInverseProblem::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The indices of the equilibrium species
    Indices ies;

    /// The indices of the equilibrium elements
    Indices iee;

    /// The number of species and elements in the system
    unsigned N, E;

    /// The number of species and elements in the equilibrium partition
    unsigned Ne, Ee;

    /// The given values of the species activities
    std::vector<double> a_values;

    /// The given values of the species amounts (in units of mol)
    std::vector<double> n_values;

    /// The indices of the species whose activities are given
    Indices a_indices;

    /// The indices of the species whose amounts are given
    Indices n_indices;

    /// The initial amounts of the elements in the equilibrium partition (in units of mol)
    Vector b0;

    /// The names of the titrants
    std::vector<std::string> titrants;

    /// The formulas of the titrants
    std::vector<std::map<std::string, double>> formulas;

    /// The pairs of mutually exclusive titrants
    std::vector<std::pair<std::string, std::string>> exclusive;

    /// The formula matrix of the titrants
    Matrix C;

    /// Construct an Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialize the number of species and elements in the system
        N = system.numSpecies();
        E = system.numElements();

        // Set the default partition as all species are in equilibrium
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical system
        partition = partition_;

        // Initialize the number of species and elements in the equilibrium partition
        Ne = partition.numEquilibriumSpecies();
        Ee = partition.numEquilibriumElements();

        // Initialize the indices of the equilibrium species and elements
        ies = partition.indicesEquilibriumSpecies();
        iee = partition.indicesEquilibriumElements();
    }

    /// Set the partition of the chemical system using a formatted string
    auto setPartition(std::string partition) -> void
    {
        setPartition(Partition(system, partition));
    }

    /// Add an activity constraint to the inverse equilibrium problem.
    auto addActivityConstraint(std::string species, double value) -> void
    {
        const Index ispecies = system.indexSpeciesWithError(species);
        a_values.push_back(value);
        a_indices.push_back(ispecies);
    }

    /// Add an amount constraint to the inverse equilibrium problem.
    auto addAmountConstraint(std::string species, double value) -> void
    {
        const Index ispecies = system.indexSpeciesWithError(species);
        n_values.push_back(value);
        n_indices.push_back(ispecies);
    }

    /// Set the known molar amounts of the elements before unknown amounts of titrants have been added.
    auto setInitialElementAmounts(const Vector& b0) -> void
    {
        this->b0 = b0;
    }

    /// Add a titrant to the inverse equilibrium problem.
    auto addTitrant(std::string titrant, std::map<std::string, double> formula) -> void
    {
        // Assert that the new titrant has not been added before
        Assert(!contained(titrant, titrants), "Could not add the titrant " + titrant +
            " to the inverse problem.", "This titrant has been added before.");

        // Update the list of titrants
        titrants.push_back(titrant);

        // Update the list of formulas
        formulas.push_back(formula);

        // Update the formula matrix of the formulas
        C.conservativeResize(Ee, formulas.size());

        // The index of the last column in the formula matrix
        Index ilast = formulas.size() - 1;

        // Set the last create column to zero
        C.col(ilast).fill(0.0);

        // Loop over all (element, stoichiometry) pairs of titrant formula
        for(auto pair : formula)
        {
            // Get the local index of current element in the equilibrium partition
            const Index ielement = partition.indexEquilibriumElement(pair.first);

            // Set the entry of the formula matrix
            C(ielement, ilast) = pair.second;
        }
    }

    /// Add two titrants that are mutually exclusive.
    auto setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void
    {
        // Common error message
        auto errormsg = "Could not set `" + titrant1 + "` and " + titrant2 + "` "
            "as mutually exclusive titrants.";

        // Check if the two titrants are different
        Assert(titrant1 != titrant2, errormsg, "They must have different identifiers.");

        // Check if the two titrants have been added before
        Assert(contained(titrant1, titrants) && contained(titrant2, titrants), errormsg,
            "At least one of them have not been added before.");

        // Update the list of mutually exclusive titrants
        exclusive.push_back({titrant1, titrant2});
    }

    /// Return the formula matrix of the formulas.
    auto formulaMatrixTitrants() const -> Matrix
    {
        return C;
    }

    /// Return the initial amounts of elements before unknown amounts of titrants are added.
    auto initialElementAmounts() const -> Vector
    {
        return b0;
    }

    /// Return the indices of the species with given activity values.
    auto indicesSpeciesActivityConstraints() const -> Indices
    {
        return a_indices;
    }

    /// Return the indices of the species with given amount values.
    auto indicesSpeciesAmountConstraints() const -> Indices
    {
        return n_indices;
    }

    /// Return the values of the activities of the species with given activities.
    auto valuesSpeciesActivityConstraints() const -> Vector
    {
        return Vector::Map(a_values.data(), a_values.size());
    }

    /// Return the values of the amounts of the species with given amounts.
    auto valuesSpeciesAmountConstraints() const -> Vector
    {
        return Vector::Map(n_values.data(), n_values.size());
    }
};

EquilibriumInverseProblem::EquilibriumInverseProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumInverseProblem::EquilibriumInverseProblem(const EquilibriumInverseProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumInverseProblem::~EquilibriumInverseProblem()
{}

auto EquilibriumInverseProblem::operator=(EquilibriumInverseProblem other) -> EquilibriumInverseProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumInverseProblem::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumInverseProblem::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumInverseProblem::addActivityConstraint(std::string species, double value) -> void
{
    pimpl->addActivityConstraint(species, value);
}

auto EquilibriumInverseProblem::addAmountConstraint(std::string species, double value) -> void
{
    pimpl->addAmountConstraint(species, value);
}

auto EquilibriumInverseProblem::setInitialElementAmounts(const Vector& b0) -> void
{
    pimpl->setInitialElementAmounts(b0);
}

auto EquilibriumInverseProblem::addTitrant(std::string titrant, std::map<std::string, double> formula) -> void
{
    pimpl->addTitrant(titrant, formula);
}

auto EquilibriumInverseProblem::setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void
{
    pimpl->setAsMutuallyExclusive(titrant1, titrant2);
}

auto EquilibriumInverseProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumInverseProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

auto EquilibriumInverseProblem::formulaMatrixTitrants() const -> Matrix
{
    return pimpl->formulaMatrixTitrants();
}

auto EquilibriumInverseProblem::initialElementAmounts() const -> Vector
{
    return pimpl->initialElementAmounts();
}

auto EquilibriumInverseProblem::indicesSpeciesActivityConstraints() const -> Indices
{
    return pimpl->indicesSpeciesActivityConstraints();
}

auto EquilibriumInverseProblem::indicesSpeciesAmountConstraints() const -> Indices
{
    return pimpl->indicesSpeciesAmountConstraints();
}

auto EquilibriumInverseProblem::valuesSpeciesActivityConstraints() const -> Vector
{
    return pimpl->valuesSpeciesActivityConstraints();
}

auto EquilibriumInverseProblem::valuesSpeciesAmountConstraints() const -> Vector
{
    return pimpl->valuesSpeciesAmountConstraints();
}

} // namespace Reaktoro
