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
#include <Reaktoro/Common/ElementUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {
namespace {

/// A type used to define the residual of an equilibrium constraint and its partial derivatives.
struct ResidualEquilibriumConstraint
{
    /// Construct a default EquilibriumConstraintResult instance
    ResidualEquilibriumConstraint() {};

    /// Construct an EquilibriumConstraintResult instance from a ChemicalScalar instance
    ResidualEquilibriumConstraint(const ChemicalScalar& scalar)
    : val(scalar.val), ddn(scalar.ddn) {}

    /// The residual value of the equilibrium constraint.
    double val;

    /// The partial derivatives of the residual w.r.t. titrant amounts x.
    Vector ddx;

    /// The partial derivatives of the residual w.r.t. species amounts n.
    Vector ddn;

    /// Assign a ChemicalScalar instance to this.
    auto operator=(const ChemicalScalar& scalar) -> ResidualEquilibriumConstraint&
    {
        val = scalar.val;
        ddn = scalar.ddn;
        return *this;
    }
};

/// A type used to define the functional signature of an equilibrium constraint.
using EquilibriumConstraint =
    std::function<ResidualEquilibriumConstraint
        (const Vector&, const ChemicalState&)>;

/// Return a Titrant instance constructed from a single component titrant or from a multi-component titrant.
auto parseTitrant(std::string titrant) -> Titrant
{
    // First, split on ( and ) to check if titrant has solution format like (1:kg:H2O)(1:mol:NaCl)
    const auto words = split(titrant, "()");

    // Check if size is 1, so that titrant is in a compound format like H2O
    if(words.size() == 1)
        return Titrant(titrant, Reaktoro::elements(titrant));

    // The elemntal formula of the titrant in solution format
    std::map<std::string, double> formula;

    // Otherwise, process each word
    for(auto word : words)
    {
        // Split on `:` to get three words, like 1:kg:H2O results in ("1", "kg", "H2O")
        const auto triplet = split(word, ":");

        // Define some auxiliary references
        const auto number = tofloat(triplet[0]);
        const auto units = triplet[1];
        const auto name = triplet[2];

        // Get the elements that compose the current titrant component
        const auto elements = Reaktoro::elements(name);

        // Get the molar mass of the current titrant component
        const auto molar_mass = Reaktoro::molarMass(elements);

        // The conversion factor from given units to mol
        double factor = 1;

        // Check if the type of units, either mass or mol
        if(units::convertible(units, "kg"))
            factor = units::convert(1.0, units, "kg")/molar_mass; // convert from mass units to kg, then from kg to mol
        else if(units::convertible(units, "mol"))
            factor = units::convert(1.0, units, "mol"); // convert from molar units to mol
        else RuntimeError("Could not create the titrant with name `" + titrant + "`.",
            "The units of each component in a multi-component titrant must be convertible to mol or kg.");

        // Determine the molar coefficient of the titrant (e.g., 1 kg H2O ~ 55.508 mol H2O)
        const auto coeff = factor * number;

        // Update the elemental formula of the titrant in solution format
        for(auto pair : elements)
            if(formula.count(pair.first))
                formula[pair.first] += coeff * pair.second;
            else
                formula.insert({pair.first, coeff * pair.second});
    }

    return Titrant(titrant, formula);
}

} // namespace

Titrant::Titrant()
: _molar_mass(0.0)
{}

Titrant::Titrant(const Species& species)
: _name(species.name()), _molar_mass(species.molarMass())
{
    for(auto x : species.elements())
        _formula.insert({x.first.name(), x.second});
}

Titrant::Titrant(std::string titrant)
: Titrant(parseTitrant(titrant))
{}

Titrant::Titrant(std::string name, const std::map<std::string, double>& formula)
: _name(name), _formula(formula), _molar_mass(Reaktoro::molarMass(formula))
{}

auto Titrant::name() const -> std::string
{
    return _name;
}

auto Titrant::formula() const -> std::map<std::string, double>
{
    return _formula;
}

auto Titrant::formula(const ChemicalSystem& system) const -> Vector
{
    Vector c = zeros(system.numElements());
    for(auto pair : _formula)
        c[system.indexElement(pair.first)] = pair.second;
    return c;
}

auto Titrant::molarMass() const -> double
{
    return _molar_mass;
}

auto operator<(const Titrant& lhs, const Titrant& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Titrant& lhs, const Titrant& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

struct EquilibriumInverseProblem::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The equilibrium constraint functions
    std::vector<EquilibriumConstraint> constraints;

    /// The initial amounts of the elements in the equilibrium partition (in units of mol)
    Vector b0;

    /// The initial guess of the titrants (in units of mol)
    Vector x0;

    /// The names of the titrants
    std::vector<Titrant> titrants;

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
    auto addSpeciesActivityConstraint(std::string species, double value) -> void
    {
        // The index of the species
        const Index ispecies = system.indexSpeciesWithError(species);

        // The ln of the given activity value
        const double ln_val = std::log(value);

        // Auxiliary chemical scalar to avoid memory reallocation
        ChemicalScalar ln_ai;

        // Define the activity constraint function
        EquilibriumConstraint f = [=](const Vector& x, const ChemicalState& state) mutable
        {
            ln_ai = state.properties().lnActivities()[ispecies];
            return ln_ai - ln_val;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Add an amount constraint to the inverse equilibrium problem.
    auto addSpeciesAmountConstraint(std::string species, double value) -> void
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
        EquilibriumConstraint f = [=](const Vector& x, const ChemicalState& state) mutable
        {
            ni.val = state.speciesAmount(ispecies);
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
        EquilibriumConstraint f = [=](const Vector& x, const ChemicalState& state) mutable
        {
            np = state.properties().phaseMoles()[iphase];
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
        EquilibriumConstraint f = [=](const Vector& x, const ChemicalState& state) mutable
        {
            Vp = state.properties().phaseVolumes()[iphase];
            return Vp - value;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Set the known molar amounts of the elements before unknown amounts of titrants have been added.
    auto setElementInitialAmounts(const Vector& binit) -> void
    {
        b0 = binit;
    }

    /// Set the initial guess of a titrant.
    auto setTitrantInitialAmount(std::string titrant, double amount) -> void
    {
        // Get the index of the titrant
        const Index ititrant = index(titrant, titrants);

        // Assert this titrant has been added already
        Assert(ititrant < titrants.size(), "Could not set the initial guess "
            "of titrant `" + titrant + "`.", "This titrant has not been added yet.");

        // Set the initial guess of the titrant
        x0[ititrant] = amount;
    }

    /// Add a titrant to the inverse equilibrium problem.
    auto addTitrant(const Titrant& titrant) -> void
    {
        // Assert that the new titrant has not been added before
        Assert(!contained(titrant, titrants), "Could not add the titrant " + titrant.name() +
            " to the inverse problem.", "This titrant has been added before.");

        // Update the list of titrants
        titrants.push_back(titrant);

        // The number of elements in the system
        const Index num_elements = system.numElements();

        // The index of the last titrant
        const Index ilast = titrants.size() - 1;

        // Update the formula matrix of the formulas
        formula_matrix_titrants.conservativeResize(num_elements, titrants.size());

        // Set the last created column to zero
        formula_matrix_titrants.col(ilast) = titrant.formula(system);

        // Resize the vector of titrant initial guess
        x0.conservativeResize(titrants.size());

        // Initialize the initial guess of the new titrant to zero
        x0[ilast] = 0.0;
    }

    /// Add a titrant to the inverse equilibrium problem using either a species name or a compound.
    auto addTitrant(std::string titrant) -> void
    {
        const Index ispecies = system.indexSpecies(titrant);
        if(ispecies < system.numSpecies())
            addTitrant(Titrant(system.species(ispecies)));
        else addTitrant(Titrant(titrant));
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

        // Get the indices of the titrants
        const Index i1 = index(titrant1, titrants);
        const Index i2 = index(titrant2, titrants);

        // The smoothing parameter for the mutually exclusive constraint function
        const double tau = 1e-20;

        // Define the mutually exclusive constraint function
        ResidualEquilibriumConstraint res;
        EquilibriumConstraint f = [=](const Vector& x, const ChemicalState& state) mutable
        {
            const Index Nt = x.rows();

            res.ddx.resize(Nt);

            const double x1 = x[i1];
            const double x2 = x[i2];

            res.val = x1*x2 - tau;
            res.ddx = x1 * unit(Nt, i2) + x2 * unit(Nt, i1);

            return res;
        };

        // Update the list of constraint functions
        constraints.push_back(f);
    }

    /// Return the formula matrix of the formulas.
    auto formulaMatrixTitrants() const -> Matrix
    {
        return formula_matrix_titrants;
    }

    /// Return the initial amounts of elements before unknown amounts of titrants are added.
    auto elementInitialAmounts() const -> Vector
    {
        return b0;
    }

    /// Return the initial amounts of titrants for the inverse equilibrium calculation
    auto titrantInitialAmounts() const -> Vector
    {
        return x0;
    }

    /// Return the residual of the equilibrium constraints and their partial molar derivatives.
    auto residualEquilibriumConstraints(const Vector& x, const ChemicalState& state) const -> ResidualEquilibriumConstraints
    {
        const Index num_species = system.numSpecies();
        const Index num_constraints = constraints.size();
        const Index num_titrants = titrants.size();
        ResidualEquilibriumConstraints res;
        res.val.resize(num_constraints);
        res.ddx.resize(num_constraints, num_titrants);
        res.ddn.resize(num_constraints, num_species);
        ResidualEquilibriumConstraint aux;
        for(Index i = 0; i < num_constraints; ++i)
        {
            aux = constraints[i](x, state);
            res.val[i]     = aux.val;
            res.ddx.row(i) = aux.ddx.size() ? aux.ddx : zeros(num_titrants);
            res.ddn.row(i) = aux.ddn.size() ? aux.ddn : zeros(num_species);
        }
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

auto EquilibriumInverseProblem::addSpeciesActivityConstraint(std::string species, double value) -> void
{
    pimpl->addSpeciesActivityConstraint(species, value);
}

auto EquilibriumInverseProblem::addSpeciesAmountConstraint(std::string species, double value) -> void
{
    pimpl->addSpeciesAmountConstraint(species, value);
}

auto EquilibriumInverseProblem::addPhaseAmountConstraint(std::string phase, double value) -> void
{
    pimpl->addPhaseAmountConstraint(phase, value);
}

auto EquilibriumInverseProblem::addPhaseVolumeConstraint(std::string phase, double value) -> void
{
    pimpl->addPhaseVolumeConstraint(phase, value);
}

auto EquilibriumInverseProblem::setElementInitialAmounts(const Vector& b0) -> void
{
    pimpl->setElementInitialAmounts(b0);
}

auto EquilibriumInverseProblem::setTitrantInitialAmount(std::string titrant, double amount) -> void
{
    pimpl->setTitrantInitialAmount(titrant, amount);
}

auto EquilibriumInverseProblem::addTitrant(const Titrant& titrant) -> void
{
    pimpl->addTitrant(titrant);
}

auto EquilibriumInverseProblem::addTitrant(std::string titrant) -> void
{
    pimpl->addTitrant(titrant);
}

auto EquilibriumInverseProblem::setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void
{
    pimpl->setAsMutuallyExclusive(titrant1, titrant2);
}

auto EquilibriumInverseProblem::empty() const -> bool
{
    return pimpl->constraints.empty();
}

auto EquilibriumInverseProblem::numConstraints() const -> Index
{
    return pimpl->constraints.size();
}

auto EquilibriumInverseProblem::numTitrants() const -> Index
{
    return pimpl->titrants.size();
}

auto EquilibriumInverseProblem::formulaMatrixTitrants() const -> Matrix
{
    return pimpl->formulaMatrixTitrants();
}

auto EquilibriumInverseProblem::elementInitialAmounts() const -> Vector
{
    return pimpl->elementInitialAmounts();
}

auto EquilibriumInverseProblem::titrantInitialAmounts() const -> Vector
{
    return pimpl->titrantInitialAmounts();
}

auto EquilibriumInverseProblem::residualEquilibriumConstraints(const Vector& x, const ChemicalState& state) const -> ResidualEquilibriumConstraints
{
    return pimpl->residualEquilibriumConstraints(x, state);
}

} // namespace Reaktoro
