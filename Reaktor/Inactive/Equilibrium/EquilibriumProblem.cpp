/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "EquilibriumProblem.hpp"

// C++ includes
#include <map>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Equilibrium/BalanceConstraints.hpp>
#include <Reaktor/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktor/Utils/ElementUtils.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>
#include <Reaktor/Utils/SetUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {
namespace internal {

inline auto nonAmountOrMassUnitError(std::string unit) -> void
{
    Exception exception;
    exception.error << "Cannot set the amount or mass of a species.";
    exception.reason << "The provided unit " << unit << " is not convertible to units of amount or mass (e.g., mol and g).";
    raise(exception);
}

inline auto infoAcidity(double value) -> std::string
{
    std::stringstream ss;
    ss << "pH of set to " << value;
    return ss.str();
}

inline auto infoSpeciesActivity(const std::string& species, double value) -> std::string
{
    std::stringstream ss;
    ss << "activity of species " << species << " set to " << value;
    return ss.str();
}

inline auto infoElementAmount(const std::string& element, double value) -> std::string
{
    std::stringstream ss;
    ss << "amount of element " << element << " set to " << value << " moles";
    return ss.str();
}

inline auto infoSpeciesAmount(const std::string& species, units::Amount value) -> std::string
{
    std::stringstream ss;
    ss << "amount of species " << species << " set to " << value;
    return ss.str();
}

inline auto infoSpeciesMass(const std::string& species, units::Mass value) -> std::string
{
    std::stringstream ss;
    ss << "mass of species " << species << " set to " << value;
    return ss.str();
}

inline auto infoPartialPressure(const std::string& species, units::Pressure value) -> std::string
{
    std::stringstream ss;
    ss << "partial pressure of species " << species << " set to " << value;
    return ss.str();
}

inline auto infoChargeBalance(const std::string& phase) -> std::string
{
    std::stringstream ss;
    ss << "charge balance of phase " << phase;
    return ss.str();
}

} /* namespace internal */

using namespace internal;

/**
 * Type that defines the function signature of a scalar-valued constraint function
 */
using EquilibriumConstraint = std::function<PartialScalar(const ChemicalState&)>;

class EquilibriumProblem::Impl
{
    /// The chemical system instance
    ChemicalSystem system$;

    /// The partitioning of the species in the chemical system
    Partitioning partitioning$;

    /// The mass-balance and charge-balance constraints
    BalanceConstraints balance$;

    /// The molar abundance of the elements in the equilibrium species
    std::map<std::string, double> abundance;

    /// The elements that are considered free
    std::vector<std::string> free_elements$;

    /// The non-mass-balance equilibrium constraints
    std::vector<EquilibriumConstraint> constraints$;

    /// The information of the non-mass-balance equilibrium constraints
    std::vector<std::string> info_constraints;

public:
    Impl(const ChemicalSystem& system, const Partitioning& partitioning)
    : system$(system), partitioning$(partitioning), balance$(system, partitioning)
    {}

    ~Impl() {}

    auto addCompound(std::string compound, units::Amount value) -> void
    {
        for(const auto& pair : elements(compound))
            abundance[pair.first] += pair.second * value.in(unit(mol));
    }

    auto addCompound(std::string compound, units::Mass value) -> void
    {
        units::Amount moles = value / molarMass(compound);
        addCompound(compound, moles);
    }

    auto addCompound(std::string compound, double value, std::string unit) -> void
    {
        if(units::convertible(unit, "mol"))
            addCompound(compound, units::Amount(value, unit));
        else if(units::convertible(unit, "g"))
            addCompound(compound, units::Mass(value, unit));
        else nonAmountOrMassUnitError(unit);
    }

    auto addSpecies(const Index& idx_species, units::Amount value) -> void
    {
        const auto& species = system$.species(idx_species);
        const auto& formula = species.elements();

        for(const auto& pair : formula)
            abundance[pair.first] += pair.second * value.in(unit(mol));
    }

    auto addSpecies(const Index& idx_species, units::Mass value) -> void
    {
        const auto molar_mass = system$.species(idx_species).molarMass();
        const units::Amount moles = value/molar_mass;
        addSpecies(idx_species, moles);
    }

    auto addSpecies(std::string species, units::Amount value) -> void
    {
        addSpecies(system$.idxSpeciesWithError(species), value);
    }

    auto addSpecies(std::string species, units::Mass value) -> void
    {
        addSpecies(system$.idxSpeciesWithError(species), value);
    }

    auto addSpecies(std::string species, double value, std::string unit) -> void
    {
        if(units::convertible(unit, "mol"))
            addSpecies(species, units::Amount(value, unit));
        else if(units::convertible(unit, "g"))
            addSpecies(species, units::Mass(value, unit));
        else nonAmountOrMassUnitError(unit);
    }

    auto add(const ChemicalState& state) -> void
    {
        for(const auto& idxspecies : partitioning$.idxEquilibriumSpecies())
            addSpecies(idxspecies, state.amount(idxspecies));
    }

    auto add(std::string entity, double value, std::string unit) -> void
    {
        if(system$.containsSpecies(entity)) addSpecies(entity, value, unit);
        else addCompound(entity, value, unit);
    }

    auto setAcidity(double value) -> void
    {
        constraints$.push_back(constraintAcidity(value));
        info_constraints.push_back(infoAcidity(value));
    }

    auto setChargeBalance() -> void
    {
        setChargeBalance("Aqueous");
    }

    auto setChargeBalance(std::string phase) -> void
    {
        constraints$.push_back(constraintChargeBalance(phase));
        info_constraints.push_back(infoChargeBalance(phase));
    }

    auto setSpecies(std::string species, units::Amount value) -> void
    {
        constraints$.push_back(constraintSpeciesAmount(species, value));
        info_constraints.push_back(infoSpeciesAmount(species, value));
    }

    auto setSpecies(std::string species, units::Mass value) -> void
    {
        constraints$.push_back(constraintSpeciesMass(species, value));
        info_constraints.push_back(infoSpeciesMass(species, value));
    }

    auto setActivity(std::string species, double value) -> void
    {
        constraints$.push_back(constraintActivity(species, value));
        info_constraints.push_back(infoSpeciesActivity(species, value));
    }

    auto setPartialPressure(std::string species, units::Pressure value) -> void
    {
        constraints$.push_back(constraintPartialPressure(species, value.in(unit(Pa))));
        info_constraints.push_back(infoPartialPressure(species, value));
    }

    auto setFreeElements(std::string elements) -> void
    {
        free_elements$ = split(elements, " ");
    }

    auto setFreeElements(std::vector<std::string> elements) -> void
    {
        free_elements$ = elements;
    }

    auto freeElements() const -> const std::vector<std::string>&
    {
        return free_elements$;
    }

    auto infoConstraints() -> std::vector<std::string>
    {
        return info_constraints;
    }

    auto constraints() const -> EquilibriumConstraints
    {
        if(free_elements$.empty())
        {
            Vector be = zeros(partitioning$.numEquilibriumElements());
            for(const auto& pair : abundance)
                be[partitioning$.idxEquilibriumElement(pair.first)] = pair.second;

            return balance$.constraints(be);;
        }
        else
        {
            // Create the remaining constraints that impose the number of moles of the elements
            std::vector<EquilibriumConstraint> constraints = constraints$;
            for(const auto& pair : abundance)
            {
                const auto element = pair.first;
                const auto moles = pair.second;
                if(not contained(element, free_elements$))
                    constraints.push_back(constraintElementMoles(element, moles));
            }

            // The number of constraints and equilibrium species
            const auto Nc = constraints.size();
            const auto Ne = partitioning$.numEquilibriumSpecies();

            // Define some auxiliary variables
            PartialScalar aux;
            PartialVector he;
            func(he).resize(Nc);
            grad(he).resize(Nc, Ne);

            // Define the vector-valued equilibrium constraint function
            auto fn = [=](const ChemicalState& state) mutable -> PartialVector
            {
                // Iterate over all constraints an evaluate them
                for(unsigned i = 0; i < Nc; ++i)
                {
                    aux = constraints[i](state);
                    func(he)[i] = func(aux);
                    grad(he).row(i) = grad(aux);
                }

                return he;
            };

            return EquilibriumConstraints(fn, Nc);
        }
    }

    auto constraintSpeciesAmount(const std::string& species, units::Amount value) const -> EquilibriumConstraint
    {
        const auto Ne          = partitioning$.numEquilibriumSpecies();
        const auto idx_species = partitioning$.idxEquilibriumSpecies(species);
        Vector ne;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            func(he) = ne[idx_species] - value.in(unit(mol));
            grad(he) = Vector::Unit(Ne, idx_species);
            return he;
        };

        return constraint;
    }

    auto constraintSpeciesMass(const std::string& species, units::Mass value) const -> EquilibriumConstraint
    {
        const auto idx_species = system$.idxSpeciesWithError(species);
        const auto molar_mass = system$.species(idx_species).molarMass();

        units::Amount moles = value/molar_mass;

        return constraintSpeciesAmount(species, moles);
    }

    auto constraintSpeciesMolality(const std::string& species, double molality) const -> EquilibriumConstraint
    {
        const auto Ne          = partitioning$.numEquilibriumSpecies();
        const auto idx_species = partitioning$.idxEquilibriumSpecies(species);
        const auto idx_water   = partitioning$.idxEquilibriumSpecies("H2O(l)");
        Vector ne;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            func(he) = ne[idx_species] - molality/55.508 * ne[idx_water];
            grad(he) = zeros(Ne);
            grad(he)[idx_species] = 1.0;
            grad(he)[idx_water]   = -molality/55.508;
            return he;
        };

        return constraint;
    }

    auto constraintSpeciesMolarFraction(const std::string& species, double value) const -> EquilibriumConstraint
    {
        const auto Ne                = partitioning$.numEquilibriumSpecies();
        const auto idx_species       = partitioning$.idxEquilibriumSpecies(species);
        const auto idx_phase         = system$.idxPhaseWithSpecies(species);
        const auto phase_species     = system$.phase(idx_phase).speciesNames();
        const auto idx_phase_species = partitioning$.idxEquilibriumSpecies(phase_species);
        Vector ne;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            const double nt = rows(idx_phase_species, ne).sum();

            func(he) = ne[idx_species] - value * nt;
            grad(he) = zeros(Ne);
            for(const Index& i : idx_phase_species)
                grad(he)[i] = -value;
            grad(he)[idx_species] += 1.0;

            return he;
        };

        return constraint;
    }

    auto constraintActivity(const std::string& species, double value) const -> EquilibriumConstraint
    {
        const auto idx_species = partitioning$.idxEquilibriumSpecies(species);
        PartialVector ae;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ae       = state.activities();
            func(ae) = partitioning$.equilibriumRows(func(ae));
            grad(ae) = partitioning$.equilibriumRowsCols(grad(ae));
            func(he) = func(ae)[idx_species] - value;
            grad(he) = grad(ae).row(idx_species);
            return he;
        };

        return constraint;
    }

    auto constraintElementMoles(const std::string& element, double value) const -> EquilibriumConstraint
    {
        const Matrix We = partitioning$.equilibriumFormulaMatrix(system$);
        const Vector we = We.row(system$.idxElement(element));
        Vector ne;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            func(he) = we.dot(ne) - value;
            grad(he) = we;
            return he;
        };

        return constraint;
    }

    auto constraintAcidity(double value) const -> EquilibriumConstraint
    {
        return constraintActivity("H+", std::pow(10.0, -value));
    }

    auto constraintPartialPressure(const std::string& gas, double value) const -> EquilibriumConstraint
    {
        const auto Ne                = partitioning$.numEquilibriumSpecies();
        const auto idx_gas           = partitioning$.idxEquilibriumSpecies(gas);
        const auto idx_phase         = system$.idxPhase("Gaseous");
        const auto phase_species     = system$.phase(idx_phase).speciesNames();
        const auto idx_phase_species = partitioning$.idxEquilibriumSpecies(phase_species);
        Vector ne;
        PartialScalar he;

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            const double P  = state.pressure();
            const double nt = rows(idx_phase_species, ne).sum();

            func(he) = ne[idx_gas] - value/P * nt;
            grad(he) = zeros(Ne);
            for(const Index& i : idx_phase_species)
                grad(he)[i] = -value/P;
            grad(he)[idx_gas] += 1.0;

            return he;
        };

        return constraint;
    }

    auto constraintChargeBalance(const std::string& name) const -> EquilibriumConstraint
    {
        const auto Ne                = partitioning$.numEquilibriumSpecies();
        const auto idx_phase         = system$.idxPhase(name);
        const auto phase_species     = system$.phase(idx_phase).speciesNames();
        const auto idx_phase_species = partitioning$.idxEquilibriumSpecies(phase_species);
        Vector ne;
        PartialScalar he;

        const Phase& phase = system$.phase(idx_phase);
        Vector ze = zeros(Ne);
        for(unsigned i = 0; i < idx_phase_species.size(); ++i)
            ze[idx_phase_species[i]] = phase.species(i).charge();

        EquilibriumConstraint constraint = [=](const ChemicalState& state) mutable
        {
            ne = state.equilibriumComposition(partitioning$);
            func(he) = ne.dot(ze);
            grad(he) = ze;
            return he;
        };

        return constraint;
    }
};

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: pimpl(new Impl(system, Partitioning(system)))
{}

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system, const Partitioning& partitioning)
: pimpl(new Impl(system, partitioning))
{}

EquilibriumProblem::EquilibriumProblem(const EquilibriumProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumProblem::~EquilibriumProblem()
{}

auto EquilibriumProblem::operator=(EquilibriumProblem other) -> EquilibriumProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumProblem::addCompound(std::string compound, units::Amount value) -> void
{
    pimpl->addCompound(compound, value);
}

auto EquilibriumProblem::addCompound(std::string compound, units::Mass value) -> void
{
    pimpl->addCompound(compound, value);
}

auto EquilibriumProblem::addCompound(std::string compound, double value, std::string unit) -> void
{
    pimpl->addCompound(compound, value, unit);
}

auto EquilibriumProblem::addSpecies(Index idx, units::Amount value) -> void
{
    pimpl->addSpecies(idx, value);
}

auto EquilibriumProblem::addSpecies(Index idx, units::Mass value) -> void
{
    pimpl->addSpecies(idx, value);
}

auto EquilibriumProblem::addSpecies(std::string name, units::Amount value) -> void
{
    pimpl->addSpecies(name, value);
}

auto EquilibriumProblem::addSpecies(std::string name, units::Mass value) -> void
{
    pimpl->addSpecies(name, value);
}

auto EquilibriumProblem::addSpecies(std::string species, double value, std::string unit) -> void
{
    pimpl->addSpecies(species, value, unit);
}

auto EquilibriumProblem::add(const ChemicalState& state) -> void
{
    pimpl->add(state);
}

auto EquilibriumProblem::add(std::string entity, double value, std::string unit) -> void
{
    pimpl->add(entity, value, unit);
}

auto EquilibriumProblem::setAcidity(double value) -> void
{
    pimpl->setAcidity(value);
}

auto EquilibriumProblem::setChargeBalance() -> void
{
    pimpl->setChargeBalance();
}

auto EquilibriumProblem::setChargeBalance(std::string phase) -> void
{
    pimpl->setChargeBalance(phase);
}

auto EquilibriumProblem::setSpecies(std::string species, units::Mass value) -> void
{
    pimpl->setSpecies(species, value);
}

auto EquilibriumProblem::setSpecies(std::string species, units::Amount value) -> void
{
    pimpl->setSpecies(species, value);
}

auto EquilibriumProblem::setActivity(std::string species, double value) -> void
{
    pimpl->setActivity(species, value);
}

auto EquilibriumProblem::setPartialPressure(std::string gas, units::Pressure value) -> void
{
    pimpl->setPartialPressure(gas, value);
}

auto EquilibriumProblem::setFreeElements(std::string elements) -> void
{
    pimpl->setFreeElements(elements);
}

auto EquilibriumProblem::setFreeElements(std::vector<std::string> elements) -> void
{
    pimpl->setFreeElements(elements);
}

auto EquilibriumProblem::freeElements() const -> const std::vector<std::string>&
{
    return pimpl->freeElements();
}

auto EquilibriumProblem::infoConstraints() -> std::vector<std::string>
{
    return pimpl->infoConstraints();
}

auto EquilibriumProblem::constraints() const -> EquilibriumConstraints
{
    return pimpl->constraints();
}

EquilibriumProblem::operator EquilibriumConstraints() const
{
    return pimpl->constraints();
}

} // namespace Reaktor
