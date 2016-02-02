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
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {
namespace {

/// A type used to define the functional signature of a general equilibrium constraint.
using EquilibriumConstraint = std::function<ChemicalScalar(const ChemicalProperties&)>;

} // namespace

struct EquilibriumInverseProblem::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The equilibrium constraint functions
    std::vector<EquilibriumConstraint> constraints;

    /// The initial amounts of the elements in the equilibrium partition (in units of mol)
    Vector b0;

    /// The names of the titrants
    std::vector<std::string> titrants;

    /// The formulas of the titrants
    std::vector<std::map<std::string, double>> formulas;

    /// The pairs of mutually exclusive titrants
    std::vector<std::pair<std::string, std::string>> exclusive;

    /// The formula matrix of the titrants
    Matrix formula_matrix_titrants;

    /// Construct an Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
    }

    /// Add an activity constraint to the inverse equilibrium problem.
    auto addActivityConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // The ln of the given activity value
        const double ln_val = std::log(value);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ln_ai;

        // Define the activity constraint function
        EquilibriumConstraint f = [=](const ChemicalProperties& properties) mutable
        {
            ln_ai = properties.lnActivities()[ispecies];
            return ln_ai - ln_val;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add an amount constraint to the inverse equilibrium problem.
    auto addAmountConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // The number of species in the system
        const Index num_species = system.numSpecies();

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ni(num_species);

        // Set the parial molar derivative of `ni`
        ni.ddn[ispecies] = 1.0;

        // Define the amount constraint function
        EquilibriumConstraint f = [=](const ChemicalProperties& properties) mutable
        {
            ni.val = properties.composition()[ispecies];
            return ni - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a phase amount constraint to the inverse equilibrium problem.
    auto addPhaseAmountConstraint(std::string phase, double value) -> void
    {
        // The index of the species
        const Index iphase = system.indexPhaseWithError(phase);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar np;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](const ChemicalProperties& properties) mutable
        {
            np = properties.phaseMoles()[iphase];
            return np - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add a phase volume constraint to the inverse equilibrium problem.
    auto addPhaseVolumeConstraint(std::string phase, double value) -> void
    {
        // The index of the species
        const Index iphase = system.indexPhaseWithError(phase);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar Vp;

        // Define the phase volume constraint function
        EquilibriumConstraint f = [=](const ChemicalProperties& properties) mutable
        {
            Vp = properties.phaseVolumes()[iphase];
            return Vp - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Set the known molar amounts of the elements before unknown amounts of titrants have been added.
    auto setInitialElementAmounts(const Vector& binit) -> void
    {
        b0 = binit;
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

        // The number of elements in the system
        const Index num_elements = system.numElements();

        // Update the formula matrix of the formulas
        formula_matrix_titrants.conservativeResize(num_elements, titrants.size());

        // The index of the last column in the formula matrix
        Index ilast = titrants.size() - 1;

        // Set the last created column to zero
        formula_matrix_titrants.col(ilast).fill(0.0);

        // Loop over all (element, stoichiometry) pairs of titrant formula
        for(auto pair : formula)
        {
            // Get the index of current element in the system
            const Index ielement = system.indexElement(pair.first);

            // Set the entry of the formula matrix
            formula_matrix_titrants(ielement, ilast) = pair.second;
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
        return formula_matrix_titrants;
    }

    /// Return the initial amounts of elements before unknown amounts of titrants are added.
    auto initialElementAmounts() const -> Vector
    {
        return b0;
    }

    /// Return the residual of the equilibrium constraints and their partial molar derivatives.
    auto residualConstraints(const ChemicalProperties& properties) const -> ChemicalVector
    {
        const Index num_species = system.numSpecies();
        const Index num_constraints = constraints.size();
        ChemicalVector res(num_constraints, num_species);
        for(Index i = 0; i < num_constraints; ++i)
            res[i] = constraints[i](properties);
        return res;
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

auto EquilibriumInverseProblem::addActivityConstraint(std::string species, double value) -> void
{
    pimpl->addActivityConstraint(species, value);
}

auto EquilibriumInverseProblem::addAmountConstraint(std::string species, double value) -> void
{
    pimpl->addAmountConstraint(species, value);
}

auto EquilibriumInverseProblem::addPhaseAmountConstraint(std::string phase, double value) -> void
{
    pimpl->addPhaseAmountConstraint(phase, value);
}

auto EquilibriumInverseProblem::addPhaseVolumeConstraint(std::string phase, double value) -> void
{
    pimpl->addPhaseVolumeConstraint(phase, value);
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

auto EquilibriumInverseProblem::formulaMatrixTitrants() const -> Matrix
{
    return pimpl->formulaMatrixTitrants();
}

auto EquilibriumInverseProblem::initialElementAmounts() const -> Vector
{
    return pimpl->initialElementAmounts();
}
auto EquilibriumInverseProblem::residualConstraints(const ChemicalProperties& properties) const -> ChemicalVector
{
    return pimpl->residualConstraints(properties);
}

} // namespace Reaktoro
