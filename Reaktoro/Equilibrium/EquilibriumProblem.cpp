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

#include "EquilibriumProblem.hpp"

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/ElementUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Math/MathUtils.hpp>

namespace Reaktoro {
namespace {

/// Throw a non-amount or non-mass error if units no compatible
auto errorNonAmountOrMassUnits(std::string units) -> void
{
    Exception exception;
    exception.error << "Cannot set the amount of a species or element.";
    exception.reason << "The provided units `" << units << "` is not convertible to units of amount or mass (e.g., mol and g).";
    RaiseError(exception);
}

} // namespace

struct EquilibriumProblem::Impl
{
    /// The reference to the ChemicalSystem instance
    ChemicalSystem system;

    /// The reference to the Partition instance
    Partition partition;

    /// The temperature for the equilibrium problem (in units of K)
    double T = 298.15;

    /// The pressure for the equilibrium problem (in units of Pa)
    double P = 1.0e+5;

    /// The amounts of the elements for the equilibrium problem (in units of mol)
    Vector b;

    /// The inverse equilibrium problem definition
    EquilibriumInverseProblem inverse_problem;

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), inverse_problem(system)
    {
        // Initialize the amounts of the elements
        b = zeros(system.numElements());

        // Set the partition of the chemical system as all species in equilibrium
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& part) -> void
    {
        partition = part;
    }
};

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
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

auto EquilibriumProblem::setPartition(const Partition& partition) -> EquilibriumProblem&
{
    pimpl->setPartition(partition);
    return *this;
}

auto EquilibriumProblem::setTemperature(double val) -> EquilibriumProblem&
{
    Assert(val > 0.0, "Cannot set temperature of the equilibrium problem.",
        "Given value must be positive.");
    pimpl->T = val;
    return *this;
}

auto EquilibriumProblem::setTemperature(double val, std::string units) -> EquilibriumProblem&
{
    return setTemperature(units::convert(val, units, "kelvin"));
}

auto EquilibriumProblem::setPressure(double val) -> EquilibriumProblem&
{
    Assert(val > 0.0, "Cannot set pressure of the equilibrium problem.",
        "Given value must be positive.");
    pimpl->P = val;
    return *this;
}

auto EquilibriumProblem::setPressure(double val, std::string units) -> EquilibriumProblem&
{
    return setPressure(units::convert(val, units, "pascal"));
}

auto EquilibriumProblem::setElementAmounts(const Vector& b) -> EquilibriumProblem&
{
    pimpl->b = b;
    return *this;
}

auto EquilibriumProblem::setElementAmounts(double amount) -> EquilibriumProblem&
{
    pimpl->b.fill(amount);
    return *this;
}

auto EquilibriumProblem::add(std::string name, double amount, std::string units) -> EquilibriumProblem&
{
    if(system().indexSpecies(name) < system().numSpecies())
        return addSpecies(name, amount, units);
    return addCompound(name, amount, units);
}

auto EquilibriumProblem::add(const ChemicalState& state, double factor) -> EquilibriumProblem&
{
    return addState(state, factor);
}

auto EquilibriumProblem::addCompound(std::string name, double amount, std::string units) -> EquilibriumProblem&
{
    double molar_amount = 0.0;

    if(units::convertible(units, "mol"))
    {
        molar_amount = units::convert(amount, units, "mol");
    }
    else if(units::convertible(units, "kg"))
    {
        const double mass = units::convert(amount, units, "kg");
        const double molar_mass = molarMass(name);
        molar_amount = mass / molar_mass;
    }
    else errorNonAmountOrMassUnits(units);

    for(const auto& pair : elements(name))
    {
        const auto element = pair.first;
        const auto coeffficient = pair.second;
        const auto ielement = system().indexElement(element);
        Assert(ielement < system().numElements(),
            "Cannot add the compound `" + name + "` to the equilibrium problem.",
            "This compound has element `" + element + "`, which is not present in the chemical system. "
            "Please note that this error can happen if this element is present in different valence state. "
            "In such case, check if there is a chemical species with same chemical formula, "
            "and use its name instead (e.g., instead of SiO2, use Quartz).");
        pimpl->b[ielement] += coeffficient * molar_amount;
    }

    return *this;
}

auto EquilibriumProblem::addSpecies(std::string name, double amount, std::string units) -> EquilibriumProblem&
{
    double molar_amount = 0.0;

    const Species& species = system().species(name);

    if(units::convertible(units, "mol"))
    {
        molar_amount = units::convert(amount, units, "mol");
    }
    else if(units::convertible(units, "kg"))
    {
        const double mass = units::convert(amount, units, "kg");
        const double molar_mass = species.molarMass();
        molar_amount = mass / molar_mass;
    }
    else errorNonAmountOrMassUnits(units);

    for(const auto& pair : species.elements())
    {
        const auto element = pair.first.name();
        const auto coeffficient = pair.second;
        const auto ielement = system().indexElement(element);
        pimpl->b[ielement] += coeffficient * molar_amount;
    }

    return *this;
}

auto EquilibriumProblem::addState(const ChemicalState& state, double factor) -> EquilibriumProblem&
{
    const Indices ies = pimpl->partition.indicesEquilibriumSpecies();
    pimpl->b += factor * state.elementAmountsInSpecies(ies);
    return *this;
}

auto EquilibriumProblem::setSpeciesAmount(std::string species, double value, std::string units) -> EquilibriumProblem&
{
    return setSpeciesAmount(species, value, units, species);
}

auto EquilibriumProblem::setSpeciesAmount(std::string species, double value, std::string units, std::string titrant) -> EquilibriumProblem&
{
    if(units::convertible(units, "mol"))
    {
        value = units::convert(value, units, "mol");
    }
    else if(units::convertible(units, "kg"))
    {
        const double mass = units::convert(value, units, "kg");
        const double molar_mass = pimpl->system.species(species).molarMass();
        value = mass / molar_mass;
    }
    pimpl->inverse_problem.addSpeciesAmountConstraint(species, value);
    pimpl->inverse_problem.addTitrant(titrant);
    pimpl->inverse_problem.setTitrantInitialAmount(species, value);
    return *this;
}

auto EquilibriumProblem::setSpeciesActivity(std::string species, double value) -> EquilibriumProblem&
{
    pimpl->inverse_problem.addSpeciesActivityConstraint(species, value);
    pimpl->inverse_problem.addTitrant(species);
    return *this;
}

auto EquilibriumProblem::setSpeciesActivity(std::string species, double value, std::string titrant) -> EquilibriumProblem&
{
    pimpl->inverse_problem.addSpeciesActivityConstraint(species, value);
    pimpl->inverse_problem.addTitrant(titrant);
    return *this;
}

auto EquilibriumProblem::setSpeciesActivity(std::string species, double value, std::string titrant1, std::string titrant2) -> EquilibriumProblem&
{
    pimpl->inverse_problem.addSpeciesActivityConstraint(species, value);
    pimpl->inverse_problem.addTitrant(titrant1);
    pimpl->inverse_problem.addTitrant(titrant2);
    pimpl->inverse_problem.setAsMutuallyExclusive(titrant1, titrant2);
    pimpl->inverse_problem.setTitrantInitialAmount(titrant1, 1e-6);
    pimpl->inverse_problem.setTitrantInitialAmount(titrant2, 1e-6);
    return *this;
}

auto EquilibriumProblem::setSpeciesFugacity(std::string species, double value, std::string units) -> EquilibriumProblem&
{
    value = units::convert(value, units, "bar");
    return setSpeciesActivity(species, value);
}

auto EquilibriumProblem::setSpeciesFugacity(std::string species, double value, std::string units, std::string titrant) -> EquilibriumProblem&
{
    value = units::convert(value, units, "bar");
    return setSpeciesActivity(species, value, titrant);
}

auto EquilibriumProblem::setPhaseAmount(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumProblem&
{
    value = units::convert(value, units, "mol");
    pimpl->inverse_problem.addPhaseAmountConstraint(phase, value);
    pimpl->inverse_problem.addTitrant(titrant);
    pimpl->inverse_problem.setTitrantInitialAmount(titrant, value);
    return *this;
}

auto EquilibriumProblem::setPhaseVolume(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumProblem&
{
    value = units::convert(value, units, "m3");
    pimpl->inverse_problem.addPhaseVolumeConstraint(phase, value);
    pimpl->inverse_problem.addTitrant(titrant);
    pimpl->inverse_problem.setTitrantInitialAmount(titrant, 1000 * value); // Assume 1000 mol/m3
    return *this;
}

auto EquilibriumProblem::setSumPhaseVolumes(const std::vector<std::string>& phases, double value, std::string units, std::string titrant) -> EquilibriumProblem&
{
    value = units::convert(value, units, "m3");
    pimpl->inverse_problem.addSumPhaseVolumesConstraint(phases, value);
    pimpl->inverse_problem.addTitrant(titrant);
    pimpl->inverse_problem.setTitrantInitialAmount(titrant, 1000 * value); // Assume 1000 mol/m3
    return *this;
}

auto EquilibriumProblem::pH(double value) -> EquilibriumProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return setSpeciesActivity("H+", aHplus);
}

auto EquilibriumProblem::pH(double value, std::string titrant) -> EquilibriumProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return setSpeciesActivity("H+", aHplus, titrant);
}

auto EquilibriumProblem::pH(double value, std::string titrant1, std::string titrant2) -> EquilibriumProblem&
{
    const double aHplus = std::pow(10.0, -value);
    return setSpeciesActivity("H+", aHplus, titrant1, titrant2);
}

auto EquilibriumProblem::pe(double value) -> EquilibriumProblem&
{
    RuntimeError("Could not impose the pe of the solution.",
        "This is currently not implemented.");
    return *this;
}

auto EquilibriumProblem::pe(double value, std::string reaction) -> EquilibriumProblem&
{
    RuntimeError("Could not impose the pe of the solution.",
        "This is currently not implemented.");
    return *this;
}

auto EquilibriumProblem::Eh(double value) -> EquilibriumProblem&
{
    RuntimeError("Could not impose the Eh of the solution.",
        "This is currently not implemented.");
    return *this;
}

auto EquilibriumProblem::Eh(double value, std::string reaction) -> EquilibriumProblem&
{
    RuntimeError("Could not impose the Eh of the solution.",
        "This is currently not implemented.");
    return *this;
}

auto EquilibriumProblem::isInverseProblem() const -> bool
{
    return !pimpl->inverse_problem.empty();
}

auto EquilibriumProblem::temperature() const -> double
{
    return pimpl->T;
}

auto EquilibriumProblem::pressure() const -> double
{
    return pimpl->P;
}

auto EquilibriumProblem::elementAmounts() const -> const Vector&
{
    return pimpl->b;
}

auto EquilibriumProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

EquilibriumProblem::operator EquilibriumInverseProblem() const
{
    pimpl->inverse_problem.setElementInitialAmounts(pimpl->b);
    return pimpl->inverse_problem;
}

} // namespace Reaktoro
