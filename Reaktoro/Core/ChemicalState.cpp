// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "ChemicalState.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

struct ChemicalState::Impl
{
    /// The chemical system instance
    const ChemicalSystem system;

    /// The temperature state of the chemical system (in K)
    real T = 298.15;

    /// The pressure state of the chemical system (in Pa)
    real P = 1.0e+05;

    /// The amounts of the chemical species (in mol)
    ArrayXr n;

    /// Construct a ChemicalState::Impl instance with given chemical system.
    Impl(const ChemicalSystem& system)
    : system(system), n(ArrayXr::Zero(system.species().size()))
    {}

    auto setTemperature(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive temperature value, ", val, "K, in a ChemicalState object.");
        T = val;
    }

    auto setTemperature(real val, String unit) -> void
    {
        setTemperature(units::convert(val, unit, "K"));
    }

    auto setPressure(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure value, ", val, "Pa, in a ChemicalState object.");
        P = val;
    }

    auto setPressure(real val, String unit) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure value, ", val, unit, ", in a ChemicalState object.");
        setPressure(units::convert(val, unit, "Pa"));
    }

    auto setSpeciesAmounts(real val) -> void
    {
        error(val < 0.0, "Cannot set a negative amount value, ", val, "mol, to the species in a ChemicalState object.");
        n.fill(val);
    }

    auto setSpeciesAmounts(ArrayXrConstRef values) -> void
    {
        assert(n.size() == values.size());
        assert(values.minCoeff() >= 0.0);
        n = values;
    }

    auto setSpeciesAmount(Index ispecies, real amount) -> void
    {
        assert(ispecies < system.species().size());
        assert(amount >= 0.0);
        n[ispecies] = amount;
    }

    auto setSpeciesAmount(String name, real amount) -> void
    {
        setSpeciesAmount(system.species().index(name), amount);
    }

    auto setSpeciesAmount(Index ispecies, real amount, String unit) -> void
    {
        setSpeciesAmount(ispecies, units::convert(amount, unit, "mol"));
    }

    auto setSpeciesAmount(String name, real amount, String unit) -> void
    {
        setSpeciesAmount(system.species().index(name), amount, unit);
    }

    auto setSpeciesMass(Index ispecies, real mass) -> void
    {
        assert(ispecies < system.species().size());
        assert(mass > 0.0);
        const auto amount = mass / system.species(ispecies).molarMass();
        setSpeciesAmount(ispecies, amount);
    }

    auto setSpeciesMass(String name, real mass) -> void
    {
        setSpeciesMass(system.species().index(name), mass);
    }

    auto setSpeciesMass(Index ispecies, real mass, String unit) -> void
    {
        setSpeciesMass(ispecies, units::convert(mass, unit, "kg"));
    }

    auto setSpeciesMass(String name, real mass, String unit) -> void
    {
        setSpeciesMass(system.species().index(name), mass, unit);
    }

    auto speciesAmount(Index ispecies) const -> real
    {
        assert(ispecies < system.species().size());
        return n[ispecies];
    }

    auto speciesAmount(String name) const -> real
    {
        return speciesAmount(system.species().index(name));
    }

    auto speciesAmount(Index index, String unit) const -> real
    {
        return units::convert(speciesAmount(index), "mol", unit);
    }

    auto speciesAmount(String name, String unit) const -> real
    {
        return speciesAmount(system.species().index(name), unit);
    }

    auto speciesMass(Index ispecies) const -> real
    {
        assert(ispecies < system.species().size());
        return n[ispecies] * system.species(ispecies).molarMass();
    }

    auto speciesMass(String name) const -> real
    {
        return speciesMass(system.species().index(name));
    }

    auto speciesMass(Index index, String unit) const -> real
    {
        return units::convert(speciesMass(index), "kg", unit);
    }

    auto speciesMass(String name, String unit) const -> real
    {
        return speciesMass(system.species().index(name), unit);
    }

    auto elementAmounts() const -> ArrayXr
    {
        const auto& A = system.formulaMatrix();
        return A * n.matrix();
    }

    auto phaseProps(Index iphase) const -> ChemicalPropsPhase
    {
        const auto offset = system.phases().numSpeciesUntilPhase(iphase);
        const auto length = system.phase(iphase).species().size();
        const auto np = n.segment(offset, length);
        ChemicalPropsPhase res(system.phase(iphase));
        res.update(T, P, np);
        return res;
    }

    auto props() const -> ChemicalProps
    {
        ChemicalProps res(system);
        res.update(T, P, n);
        return res;
    }
};

ChemicalState::ChemicalState(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalState::ChemicalState(const ChemicalState& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalState::~ChemicalState()
{}

auto ChemicalState::operator=(ChemicalState other) -> ChemicalState&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalState::setTemperature(real val) -> void
{
    pimpl->setTemperature(val);
}

auto ChemicalState::setTemperature(real val, String unit) -> void
{
    pimpl->setTemperature(val, unit);
}

auto ChemicalState::setPressure(real val) -> void
{
    pimpl->setPressure(val);
}

auto ChemicalState::setPressure(real val, String unit) -> void
{
    pimpl->setPressure(val, unit);
}

auto ChemicalState::setSpeciesAmounts(real val) -> void
{
    pimpl->setSpeciesAmounts(val);
}

auto ChemicalState::setSpeciesAmounts(ArrayXrConstRef n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::setSpeciesAmount(Index ispecies, real amount) -> void
{
    pimpl->setSpeciesAmount(ispecies, amount);
}

auto ChemicalState::setSpeciesAmount(Index ispecies, real amount, String unit) -> void
{
    pimpl->setSpeciesAmount(ispecies, amount, unit);
}

auto ChemicalState::setSpeciesAmount(String name, real amount) -> void
{
    pimpl->setSpeciesAmount(name, amount);
}

auto ChemicalState::setSpeciesAmount(String name, real amount, String unit) -> void
{
    pimpl->setSpeciesAmount(name, amount, unit);
}

auto ChemicalState::setSpeciesMass(Index ispecies, real mass) -> void
{
    pimpl->setSpeciesMass(ispecies, mass);
}

auto ChemicalState::setSpeciesMass(Index ispecies, real mass, String unit) -> void
{
    pimpl->setSpeciesMass(ispecies, mass, unit);
}

auto ChemicalState::setSpeciesMass(String name, real mass) -> void
{
    pimpl->setSpeciesMass(name, mass);
}

auto ChemicalState::setSpeciesMass(String name, real mass, String unit) -> void
{
    pimpl->setSpeciesMass(name, mass, unit);
}

auto ChemicalState::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalState::temperature() const -> real
{
    return pimpl->T;
}

auto ChemicalState::pressure() const -> real
{
    return pimpl->P;
}

auto ChemicalState::speciesAmounts() const -> ArrayXrConstRef
{
    return pimpl->n;
}

auto ChemicalState::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto ChemicalState::speciesAmount(Index ispecies) const -> real
{
    return pimpl->speciesAmount(ispecies);
}

auto ChemicalState::speciesAmount(Index ispecies, String unit) const -> real
{
    return pimpl->speciesAmount(ispecies, unit);
}

auto ChemicalState::speciesAmount(String name) const -> real
{
    return pimpl->speciesAmount(name);
}

auto ChemicalState::speciesAmount(String name, String unit) const -> real
{
    return pimpl->speciesAmount(name, unit);
}

auto ChemicalState::speciesMass(Index ispecies) const -> real
{
    return pimpl->speciesMass(ispecies);
}

auto ChemicalState::speciesMass(Index ispecies, String unit) const -> real
{
    return pimpl->speciesMass(ispecies, unit);
}

auto ChemicalState::speciesMass(String name) const -> real
{
    return pimpl->speciesMass(name);
}

auto ChemicalState::speciesMass(String name, String unit) const -> real
{
    return pimpl->speciesMass(name, unit);
}

auto ChemicalState::phaseProps(Index iphase) const -> ChemicalPropsPhase
{
    return pimpl->phaseProps(iphase);
}

auto ChemicalState::props() const -> ChemicalProps
{
    return pimpl->props();
}

} // namespace Reaktoro
