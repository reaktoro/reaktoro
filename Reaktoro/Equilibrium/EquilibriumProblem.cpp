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

    /// The collected amounts of the elements (in units of mol)
    Vector b;

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const Partition& partition)
    : system(partition.system()), partition(partition)
    {
        // Initialize the amounts of the elements
        b = zeros(system.numElements());
    }

    auto setTemperature(double val) -> void
    {
        Assert(val > 0.0, "Cannot set temperature of the equilibrium problem.",
            "Given value must be positive.");
        T = val;
    }

    auto setTemperature(double val, std::string units) -> void
    {
        return setTemperature(units::convert(val, units, "kelvin"));
    }

    auto setPressure(double val) -> void
    {
        Assert(val > 0.0, "Cannot set pressure of the equilibrium problem.",
            "Given value must be positive.");
        P = val;
    }

    auto setPressure(double val, std::string units) -> void
    {
        return setPressure(units::convert(val, units, "pascal"));
    }

    auto setElementAmounts(VectorConstRef b_) -> void
    {
        Assert(b.size() == b_.size(), "Could not set the initial mole amounts of the elements.",
            "Dimension mismatch between given vector of values and number of elements.");
        Assert(b.minCoeff() >= 0.0, "Could not set the initial mole amounts of the elements.",
            "Given values must be non-negative.");
        b = b_;
    }

    auto setElementAmount(Index ielement, double amount) -> void
    {
        Assert(ielement < system.numElements(), "Could not set the initial mole amount of the given element.",
            "Dimension mismatch between given vector of values and number of elements.");
        Assert(amount > 0.0, "Cannot set element amounts of the equilibrium problem.",
            "Given value must be positive.");
        b[ielement] = amount;
    }

    auto setElementAmount(std::string element, double amount) -> void
    {
        const auto ielement = system.indexElementWithError(element);
        setElementAmount(ielement, amount);
    }

    auto setElectricalCharge(double amount) -> void
    {
        return setElementAmount("Z", amount);
    }

    auto add(std::string name, double amount, std::string units) -> void
    {
        if(system.indexSpecies(name) < system.numSpecies())
            return addSpecies(name, amount, units);
        return addCompound(name, amount, units);
    }

    auto add(const ChemicalState& state) -> void
    {
        return addState(state);
    }

    auto addCompound(std::string name, double amount, std::string units) -> void
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
            const auto ielement = system.indexElement(element);
            Assert(ielement < system.numElements(),
                "Cannot add the compound `" + name + "` to the equilibrium problem.",
                "This compound has element `" + element + "`, which is not present in the chemical system. "
                "Please note that this error can happen if this element is present in different valence state. "
                "In such case, check if there is a chemical species with same chemical formula, "
                "and use its name instead (e.g., instead of SiO2, use Quartz).");
            b[ielement] += coeffficient * molar_amount;
        }

    }

    auto addSpecies(std::string name, double amount, std::string units) -> void
    {
        double molar_amount = 0.0;

        const Species& species = system.species(name);

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
            const auto ielement = system.indexElement(element);
            b[ielement] += coeffficient * molar_amount;
        }

    }

    auto addState(const ChemicalState& state) -> void
    {
        const auto& ies = partition.indicesEquilibriumSpecies();
        b += state.elementAmountsInSpecies(ies);
    }
};

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: pimpl(new Impl(Partition(system)))
{}

EquilibriumProblem::EquilibriumProblem(const Partition& partition)
: pimpl(new Impl(partition))
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

auto EquilibriumProblem::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with EquilibriumProblem::setPartition.",
        "EquilibriumProblem::setPartition is deprecated. "
        "Use constructor EquilibriumProblem(const Partition&) instead.");
}

auto EquilibriumProblem::setTemperature(double val) -> void
{
    pimpl->setTemperature(val);
}

auto EquilibriumProblem::setTemperature(double val, std::string units) -> void
{
    pimpl->setTemperature(val, units);
}

auto EquilibriumProblem::setPressure(double val) -> void
{
    pimpl->setPressure(val);
}

auto EquilibriumProblem::setPressure(double val, std::string units) -> void
{
    pimpl->setPressure(val, units);
}

auto EquilibriumProblem::setElementAmounts(VectorConstRef b) -> void
{
    pimpl->setElementAmounts(b);
}

auto EquilibriumProblem::setElementAmount(Index ielement, double amount) -> void
{
    pimpl->setElementAmount(ielement, amount);
}

auto EquilibriumProblem::setElementAmount(std::string element, double amount) -> void
{
    pimpl->setElementAmount(element, amount);
}

auto EquilibriumProblem::setElectricalCharge(double amount) -> void
{
    pimpl->setElectricalCharge(amount);
}

auto EquilibriumProblem::add(std::string name, double amount, std::string units) -> void
{
    pimpl->add(name, amount, units);
}

auto EquilibriumProblem::add(const ChemicalState& state) -> void
{
    pimpl->add(state);
}

auto EquilibriumProblem::addCompound(std::string name, double amount, std::string units) -> void
{
    pimpl->addCompound(name, amount, units);
}

auto EquilibriumProblem::addSpecies(std::string name, double amount, std::string units) -> void
{
    pimpl->addSpecies(name, amount, units);
}

auto EquilibriumProblem::addState(const ChemicalState& state) -> void
{
    pimpl->addState(state);
}

auto EquilibriumProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

auto EquilibriumProblem::temperature() const -> double
{
    return pimpl->T;
}

auto EquilibriumProblem::pressure() const -> double
{
    return pimpl->P;
}

auto EquilibriumProblem::elementAmounts() const -> VectorConstRef
{
    return pimpl->b;
}

} // namespace Reaktoro
