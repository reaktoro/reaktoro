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

// Optima includes
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

//=================================================================================================
//
// ChemicalState
//
//=================================================================================================

struct ChemicalState::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The properties related to an equilibrium state.
    Equilibrium equilibrium;

    /// The temperature state of the chemical system (in K)
    real T = 298.15;

    /// The pressure state of the chemical system (in Pa)
    real P = 1.0e+05;

    /// The amounts of the chemical species (in mol)
    ArrayXr n;

    /// Construct a ChemicalState::Impl instance with given chemical system.
    Impl(const ChemicalSystem& system)
    : system(system), equilibrium(system),
      n(ArrayXr::Zero(system.species().size()))
    {}

    auto setTemperature(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive temperature "
            "value, ", val, "K, in a ChemicalState object.");
        T = val;
    }

    auto setTemperature(real val, String unit) -> void
    {
        setTemperature(units::convert(val, unit, "K"));
    }

    auto setPressure(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure "
            "value, ", val, "Pa, in a ChemicalState object.");
        P = val;
    }

    auto setPressure(real val, String unit) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure "
            "value, ", val, unit, ", in a ChemicalState object.");
        setPressure(units::convert(val, unit, "Pa"));
    }

    auto setSpeciesAmounts(real val) -> void
    {
        error(val < 0.0, "Cannot set a negative species "
            "amount, ", val, " mol, in a ChemicalState object.");
        n.fill(val);
    }

    auto setSpeciesAmounts(ArrayXrConstRef values) -> void
    {
        assert(n.size() == values.size());
        assert(values.minCoeff() >= 0.0);
        n = values;
    }

    auto setSpeciesAmounts(ArrayXdConstRef values) -> void
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

auto ChemicalState::setSpeciesAmounts(ArrayXdConstRef n) -> void
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

auto ChemicalState::equilibrium() const -> const Equilibrium&
{
    return pimpl->equilibrium;
}

auto ChemicalState::equilibrium() -> Equilibrium&
{
    return pimpl->equilibrium;
}

//=================================================================================================
//
// ChemicalState::Equilibrium
//
//=================================================================================================

struct ChemicalState::Equilibrium::Impl
{
    /// The number of species in the chemical system.
    const Index Nn;

    /// The number of elements in the chemical system.
    const Index Ne;

    /// The Optima::State object used for warm start Optima optimization calculations.
    Optima::State optstate;

    /// The indices of the species partitioned as (primary, secondary).
    ArrayXl ips;

    /// The number of primary species among the species.
    Index kp = 0;

    /// The indices of elements whose amounts should be positive, but given amount was less or equal to zero.
    ArrayXl isue;

    /// The indices of species that contain one or more strictly unstable elements.
    ArrayXl isus;

    /// Construct a default ChemicalState::Equilibrium::Impl instance
    Impl(const ChemicalSystem& system)
    : Nn(system.species().size()), Ne(system.elements().size() + 1)
    {}
};

ChemicalState::Equilibrium::Equilibrium(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalState::Equilibrium::Equilibrium(const ChemicalState::Equilibrium& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalState::Equilibrium::~Equilibrium()
{}

auto ChemicalState::Equilibrium::operator=(ChemicalState::Equilibrium other) -> Equilibrium&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalState::Equilibrium::setOptimaState(const Optima::State& state) -> void
{
    pimpl->optstate = state;
}

auto ChemicalState::Equilibrium::optimaState() const -> const Optima::State&
{
    return pimpl->optstate;
}

auto ChemicalState::Equilibrium::setIndicesPrimarySecondarySpecies(ArrayXlConstRef ips, Index kp) -> void
{
    pimpl->ips = ips;
    pimpl->kp = kp;
}

auto ChemicalState::Equilibrium::setIndicesStrictlyUnstableElements(ArrayXlConstRef isue) -> void
{
    pimpl->isue = isue;
}

auto ChemicalState::Equilibrium::setIndicesStrictlyUnstableSpecies(ArrayXlConstRef isus) -> void
{
    pimpl->isus = isus;
}

auto ChemicalState::Equilibrium::numPrimarySpecies() const -> Index
{
    return pimpl->kp;
}

auto ChemicalState::Equilibrium::numSecondarySpecies() const -> Index
{
    return pimpl->ips.size() - pimpl->kp;
}

auto ChemicalState::Equilibrium::indicesPrimarySpecies() const -> ArrayXlConstRef
{
    return pimpl->ips.head(numPrimarySpecies());
}

auto ChemicalState::Equilibrium::indicesSecondarySpecies() const -> ArrayXlConstRef
{
    return pimpl->ips.tail(numSecondarySpecies());
}

auto ChemicalState::Equilibrium::indicesStrictlyUnstableElements() const -> ArrayXlConstRef
{
    return pimpl->isue;
}

auto ChemicalState::Equilibrium::indicesStrictlyUnstableSpecies() const -> ArrayXlConstRef
{
    return pimpl->isus;
}

auto ChemicalState::Equilibrium::elementChemicalPotentials() const -> ArrayXdConstRef
{
    if(pimpl->optstate.ye.size())
        return pimpl->optstate.ye.head(pimpl->Ne);
    else return pimpl->optstate.ye;
}

auto ChemicalState::Equilibrium::speciesStabilities() const -> ArrayXdConstRef
{
    if(pimpl->optstate.s.size())
        return pimpl->optstate.s.head(pimpl->Nn);
    else return pimpl->optstate.s;
}

auto ChemicalState::Equilibrium::explicitTitrantAmounts() const -> ArrayXdConstRef
{
    return pimpl->optstate.p;
}

auto ChemicalState::Equilibrium::implicitTitrantAmounts() const -> ArrayXdConstRef
{
    const auto Nx = pimpl->optstate.x.size();
    const auto Nn = pimpl->Nn;
    const auto Nq = Nx - Nn;
    if(pimpl->optstate.x.size())
        return pimpl->optstate.x.tail(Nq);
    else return pimpl->optstate.x;
}

} // namespace Reaktoro
