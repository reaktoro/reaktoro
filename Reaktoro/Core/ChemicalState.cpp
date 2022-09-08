// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Optima includes
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Utils.hpp>

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

    /// The chemical properties of the system associated to this chemical state.
    ChemicalProps props;

    /// The temperature state of the chemical system (in K)
    real T = 298.15;

    /// The pressure state of the chemical system (in Pa)
    real P = 1.0e+05;

    /// The amounts of the chemical species (in mol)
    ArrayXr n;

    /// The surface areas for the existing interphase surfaces.
    ArrayXr surface_areas;

    /// Construct a ChemicalState::Impl instance with given chemical system.
    Impl(const ChemicalSystem& system)
    : system(system), equilibrium(system), props(system)
    {
        n.setConstant(system.species().size(), 1e-16); // set small positive value for initial species amounts
        surface_areas.setZero(system.reactingPhaseInterfaces().size());
    }

    auto temperature(real val) -> void
    {
        errorif(val <= 0.0, "Expecting a non-positive temperature value, but got ", val, " K.");
        T = val;
    }

    auto temperature(real val, Chars unit) -> void
    {
        temperature(units::convert(val, unit, "K"));
    }

    auto pressure(real val) -> void
    {
        errorif(val <= 0.0, "Expecting a non-positive pressure value, but got ", val, " Pa.");
        P = val;
    }

    auto pressure(real val, Chars unit) -> void
    {
        errorif(val <= 0.0, "Expecting a non-positive pressure value, but got ", val, " ", unit, ".");
        pressure(units::convert(val, unit, "Pa"));
    }

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING THE AMOUNT OR MASS OF SPECIES
    // --------------------------------------------------------------------------------------------

    auto setSpeciesAmounts(real val) -> void
    {
        errorif(val < 0.0, "It is not possible to set a negative species amount.");
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

    auto setSpeciesAmount(const StringOrIndex& species, real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        n[ispecies] = units::convert(amount, unit, "mol");
    }

    auto setSpeciesMass(const StringOrIndex& species, real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        n[ispecies] = units::convert(mass, unit, "kg") / system.species(ispecies).molarMass();
    }

    auto set(const StringOrIndex& species, real value, Chars unit) -> void
    {
        errorif(value < 0.0, "Expecting a non-negative amount/mass value, but got ", value, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        const auto numspecies = system.species().size();
        errorif(ispecies >= numspecies, "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] = amount;
    }

    auto add(const StringOrIndex& species, real value, Chars unit) -> void
    {
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        const auto numspecies = system.species().size();
        errorif(ispecies >= numspecies, "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] += amount;
        errorif(n[ispecies] < 0.0, "It is not possible to add a negative species amount (", value, " ", unit, ") that produces a negative amount for the species.");
    }

    // --------------------------------------------------------------------------------------------
    // METHODS FOR GETTING THE AMOUNT OR MASS OF SPECIES, ELEMENTS, AND CHARGE
    // --------------------------------------------------------------------------------------------

    auto speciesAmountsInPhase(const StringOrIndex& phase) const -> ArrayXrConstRef
    {
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        return n.segment(start, size);
    }

    auto speciesAmount(const StringOrIndex& species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        return n[ispecies];
    }

    auto speciesMass(const StringOrIndex& species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species), "`.");
        return n[ispecies] * system.species(ispecies).molarMass();
    }

    auto componentAmounts() const -> ArrayXr
    {
        const auto& A = system.formulaMatrix();
        return A * n.matrix();
    }

    auto elementAmounts() const -> ArrayXr
    {
        const auto& Ae = system.formulaMatrixElements();
        return Ae * n.matrix();
    }

    auto charge() const -> real
    {
        const auto& Az = system.formulaMatrixCharge();
        return (Az * n.matrix())[0];
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE AMOUNTS OF SPECIES IN THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleSpeciesAmounts(double scalar) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        n *= scalar;
    }

    auto scaleSpeciesAmounts(double scalar, const Indices& indices) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        n(indices) *= scalar;
    }

    auto scaleSpeciesAmountsInPhase(const StringOrIndex& phase, double scalar) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        n.segment(start, size) *= scalar;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE VOLUME OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto current_volume = props.volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseVolume(const StringOrIndex& phase, real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        props.update(T, P, n);
        const auto current_volume = props.phaseProps(iphase).volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_volume =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        const auto& factor = current_fluid_volume > 0.0 ? volume / current_fluid_volume : real(0.0);
        const auto& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        n(ifluidspecies) *= factor;
    }

    auto scaleSolidVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_volume =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        const auto& factor = current_solid_volume > 0.0 ? volume / current_solid_volume : real(0.0);
        const auto& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        n(isolidspecies) *= factor;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE MASS OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto current_mass = props.mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseMass(const StringOrIndex& phase, real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        props.update(T, P, n);
        const auto current_mass = props.phaseProps(iphase).mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_mass =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        const auto& factor = current_fluid_mass > 0.0 ? mass / current_fluid_mass : real(0.0);
        const auto& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        n(ifluidspecies) *= factor;
    }

    auto scaleSolidMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_mass =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        const auto& factor = current_solid_mass > 0.0 ? mass / current_solid_mass : real(0.0);
        const auto& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        n(isolidspecies) *= factor;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING/GETTING SURFACE AREAS BETWEEN PHASES
    // --------------------------------------------------------------------------------------------

    auto setSurfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2, real value, Chars unit) -> void
    {
        errorif(value < 0.0, "Expecting a non-negative surface area value, but got ", value, " ", unit);
        value = units::convert(value, unit, "m2");
        const auto isurface = surfaceIndex(phase1, phase2);
        const auto numsurfaces = system.reactingPhaseInterfaces().size();
        errorif(isurface >= numsurfaces, "Cannot set surface area for the interface between phases `", detail::stringfy(phase1), "` and `", detail::stringfy(phase2), "` because these two phases are not reacting kinetically (i.e., there are no heteroneous reactions in the chemical system in which these two phases are present).");
        surface_areas[isurface] = value;
    }

    auto setSurfaceArea(Index isurface, real value, Chars unit) -> void
    {
        const auto numsurfaces = system.reactingPhaseInterfaces().size();
        errorif(value < 0.0, "Expecting a non-negative surface area value, but got ", value, " ", unit);
        errorif(isurface >= numsurfaces, "The given surface index,", isurface, ", is out of bounds. There are only ", numsurfaces, " reacting phase interfaces in the chemical system, automatically determined from provided heterogeneous reactions.");
        value = units::convert(value, unit, "m2");
        surface_areas[isurface] = value;
    }

    auto surfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2) const -> real
    {
        const auto numsurfaces = system.reactingPhaseInterfaces().size();
        const auto isurface = surfaceIndex(phase1, phase2);
        errorif(isurface >= numsurfaces, "Cannot set surface area for the interface between phases `", detail::stringfy(phase1), "` and `", detail::stringfy(phase2), "` because these two phases are not reacting kinetically (i.e., there are no heteroneous reactions in the chemical system in which these two phases are present).");
        return surface_areas[isurface];
    }

    auto surfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2, real value, Chars unit) -> void
    {
        setSurfaceArea(phase1, phase2, value, unit);
    }

    auto surfaceArea(Index isurface) const -> real
    {
        const auto numsurfaces = system.reactingPhaseInterfaces().size();
        errorif(isurface >= numsurfaces, "The given surface index,", isurface, ", is out of bounds. There are only ", numsurfaces, " reacting phase interfaces in the chemical system, automatically determined from provided heterogeneous reactions.");
        return surface_areas[isurface];
    }

    auto surfaceIndex(const StringOrIndex& phase1, const StringOrIndex& phase2) const -> Index
    {
        const auto iphase1 = detail::resolvePhaseIndex(system, phase1);
        const auto iphase2 = detail::resolvePhaseIndex(system, phase2);
        const auto numphases = system.phases().size();
        errorif(iphase1 >= numphases, "Could not find a phase in the system with index or name `", detail::stringfy(phase1), "`.");
        errorif(iphase2 >= numphases, "Could not find a phase in the system with index or name `", detail::stringfy(phase2), "`.");
        const auto pair = iphase1 < iphase2 ? Pair<Index,Index>{iphase1, iphase2} : Pair<Index,Index>{iphase2, iphase1};
        const auto isurface = index(system.reactingPhaseInterfaces(), pair);
        return isurface;
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

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING/GETTING TEMPERATURE
// --------------------------------------------------------------------------------------------

auto ChemicalState::setTemperature(real value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::setTemperature(real value, Chars unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::temperature(real value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::temperature(real value, Chars unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::temperature() const -> real
{
    return pimpl->T;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING/GETTING PRESSURE
// --------------------------------------------------------------------------------------------

auto ChemicalState::setPressure(real value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::setPressure(real value, Chars unit) -> void
{
    pimpl->pressure(value, unit);
}

auto ChemicalState::pressure(real value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::pressure(real value, Chars unit) -> void
{
    pimpl->pressure(value, unit);
}

auto ChemicalState::pressure() const -> real
{
    return pimpl->P;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING THE AMOUNT OR MASS OF SPECIES
// --------------------------------------------------------------------------------------------

auto ChemicalState::setSpeciesAmounts(real value) -> void
{
    pimpl->setSpeciesAmounts(value);
}

auto ChemicalState::setSpeciesAmounts(ArrayXrConstRef n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::setSpeciesAmounts(ArrayXdConstRef n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::set(const StringOrIndex& species, real value, Chars unit) -> void
{
    pimpl->set(species, value, unit);
}

auto ChemicalState::add(const StringOrIndex& species, real value, Chars unit) -> void
{
    pimpl->add(species, value, unit);
}

auto ChemicalState::setSpeciesAmount(const StringOrIndex& species, real amount, Chars unit) -> void
{
    pimpl->setSpeciesAmount(species, amount, unit);
}

auto ChemicalState::setSpeciesMass(const StringOrIndex& species, real mass, Chars unit) -> void
{
    pimpl->setSpeciesMass(species, mass, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS FOR GETTING THE AMOUNT OR MASS OF SPECIES, ELEMENTS, AND CHARGE
// --------------------------------------------------------------------------------------------

auto ChemicalState::speciesAmounts() const -> ArrayXrConstRef
{
    return pimpl->n;
}

auto ChemicalState::speciesAmountsInPhase(const StringOrIndex& phase) const -> ArrayXrConstRef
{
    return pimpl->speciesAmountsInPhase(phase);
}

auto ChemicalState::speciesAmount(const StringOrIndex& species) const -> real
{
    return pimpl->speciesAmount(species);
}

auto ChemicalState::speciesMass(const StringOrIndex& species) const -> real
{
    return pimpl->speciesMass(species);
}

auto ChemicalState::componentAmounts() const -> ArrayXr
{
    return pimpl->componentAmounts();
}

auto ChemicalState::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto ChemicalState::charge() const -> real
{
    return pimpl->charge();
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE AMOUNTS OF SPECIES IN THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleSpeciesAmounts(real scalar) -> void
{
    pimpl->scaleSpeciesAmounts(scalar);
}

auto ChemicalState::scaleSpeciesAmounts(real scalar, const Indices& indices) -> void
{
    pimpl->scaleSpeciesAmounts(scalar, indices);
}

auto ChemicalState::scaleSpeciesAmountsInPhase(const StringOrIndex& phase, real scalar) -> void
{
    pimpl->scaleSpeciesAmountsInPhase(phase, scalar);
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE VOLUME OF THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleVolume(real value, Chars unit) -> void
{
    pimpl->scaleVolume(value, unit);
}

auto ChemicalState::scalePhaseVolume(const StringOrIndex& phase, real value, Chars unit) -> void
{
    pimpl->scalePhaseVolume(phase, value, unit);
}

auto ChemicalState::scaleFluidVolume(real value, Chars unit) -> void
{
    pimpl->scaleFluidVolume(value, unit);
}

auto ChemicalState::scaleSolidVolume(real value, Chars unit) -> void
{
    pimpl->scaleSolidVolume(value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE MASS OF THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleMass(real value, Chars unit) -> void
{
    pimpl->scaleMass(value, unit);
}

auto ChemicalState::scalePhaseMass(const StringOrIndex& phase, real value, Chars unit) -> void
{
    pimpl->scalePhaseMass(phase, value, unit);
}

auto ChemicalState::scaleFluidMass(real value, Chars unit) -> void
{
    pimpl->scaleFluidMass(value, unit);
}

auto ChemicalState::scaleSolidMass(real value, Chars unit) -> void
{
    pimpl->scaleSolidMass(value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING/GETTING SURFACE AREAS BETWEEN PHASES
// --------------------------------------------------------------------------------------------

auto ChemicalState::setSurfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2, real value, Chars unit) -> void
{
    pimpl->setSurfaceArea(phase1, phase2, value, unit);
}

auto ChemicalState::setSurfaceArea(Index isurface, real value, Chars unit) -> void
{
    pimpl->setSurfaceArea(isurface, value, unit);
}

auto ChemicalState::surfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2, real value, Chars unit) -> void
{
    pimpl->surfaceArea(phase1, phase2, value, unit);
}

auto ChemicalState::surfaceArea(const StringOrIndex& phase1, const StringOrIndex& phase2) const -> real
{
    return pimpl->surfaceArea(phase1, phase2);
}

auto ChemicalState::surfaceArea(Index isurface) const -> real
{
    return pimpl->surfaceArea(isurface);
}

auto ChemicalState::surfaceAreas() const -> ArrayXrConstRef
{
    return pimpl->surface_areas;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR UPDATING CHEMICAL STATE AND ITS PROPERTIES
// --------------------------------------------------------------------------------------------

auto ChemicalState::update(const real& T, const real& P, ArrayXrConstRef n) -> void
{
    setTemperature(T);
    setPressure(P);
    setSpeciesAmounts(n);
    props().update(T, P, n);
}

auto ChemicalState::updateIdeal(const real& T, const real& P, ArrayXrConstRef n) -> void
{
    setTemperature(T);
    setPressure(P);
    setSpeciesAmounts(n);
    props().updateIdeal(T, P, n);
}

// --------------------------------------------------------------------------------------------
// MISCELLANEOUS METHODS
// --------------------------------------------------------------------------------------------

auto ChemicalState::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalState::props() const -> const ChemicalProps&
{
    return pimpl->props;
}

auto ChemicalState::props() -> ChemicalProps&
{
    return pimpl->props;
}

auto ChemicalState::equilibrium() const -> const Equilibrium&
{
    return pimpl->equilibrium;
}

auto ChemicalState::equilibrium() -> Equilibrium&
{
    return pimpl->equilibrium;
}

auto ChemicalState::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ChemicalState::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
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

    /// The number of components in the equilibrium state.
    const Index Nb;

    /// The names of the the input variables *w* used in the equilibrium calculation.
    Strings inputs;

    /// The values of the input variables *w* used in the equilibrium calculation.
    ArrayXd w;

    /// The initial component amounts in the equilibrium calculation.
    ArrayXd b;

    /// The computed control variables *q* in the equilibrium calculation.
    ArrayXd q;

    /// The computed control variables *p* in the equilibrium calculation.
    ArrayXd p;

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
    : Nn(system.species().size()), Nb(system.elements().size() + 1)
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

auto ChemicalState::Equilibrium::setInputNames(const Strings& names) -> void
{
    pimpl->inputs = names;
}

auto ChemicalState::Equilibrium::setInputValues(VectorXdConstRef w) -> void
{
    pimpl->w = w;
}

auto ChemicalState::Equilibrium::setInitialComponentAmounts(ArrayXdConstRef b) -> void
{
    pimpl->b = b;
}

auto ChemicalState::Equilibrium::setControlVariablesP(ArrayXdConstRef p) -> void
{
    pimpl->p = p;
    pimpl->optstate.p = -p;
}

auto ChemicalState::Equilibrium::setControlVariablesQ(ArrayXdConstRef q) -> void
{
    const auto Nq = pimpl->optstate.x.size() - pimpl->Nn;
    pimpl->q = q;
    if(Nq > 0)
        pimpl->optstate.x.tail(Nq) = -q;
}

auto ChemicalState::Equilibrium::setOptimaState(const Optima::State& state) -> void
{
    pimpl->optstate = state;
    const auto Nq = state.x.size() - pimpl->Nn;
    pimpl->q = -state.x.tail(Nq); // negative because of the preferred positive coefficients in the conservation matrix
    pimpl->p = -state.p; // negative because of the preferred positive coefficients in the conservation matrix
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
        return pimpl->optstate.ye.head(pimpl->Nb);
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
    return p();
}

auto ChemicalState::Equilibrium::implicitTitrantAmounts() const -> ArrayXdConstRef
{
    return q();
}

auto ChemicalState::Equilibrium::inputNames() const -> const Strings&
{
    return pimpl->inputs;
}

auto ChemicalState::Equilibrium::inputValues() const -> VectorXdConstRef
{
    return pimpl->w;
}

auto ChemicalState::Equilibrium::initialComponentAmounts() const -> ArrayXdConstRef
{
    return pimpl->b;
}

auto ChemicalState::Equilibrium::controlVariablesP() const -> ArrayXdConstRef
{
    return pimpl->p;
}

auto ChemicalState::Equilibrium::controlVariablesQ() const -> ArrayXdConstRef
{
    return pimpl->q;
}

auto ChemicalState::Equilibrium::p() const -> ArrayXdConstRef
{
    return pimpl->p;
}

auto ChemicalState::Equilibrium::q() const -> ArrayXdConstRef
{
    return pimpl->q;
}

auto ChemicalState::Equilibrium::w() const -> VectorXdConstRef
{
    return pimpl->w;
}

auto ChemicalState::Equilibrium::b() const -> ArrayXdConstRef
{
    return pimpl->b;
}

auto ChemicalState::Equilibrium::optimaState() const -> const Optima::State&
{
    return pimpl->optstate;
}

auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&
{
    const auto& n = state.speciesAmounts();
    const auto& b = state.elementAmounts();
    const auto& phases = state.system().phases();
    const auto& species = state.system().species();
    const auto& elements = state.system().elements();
    const auto& surfaces = state.system().reactingPhaseInterfaces();
    const auto& surface_areas = state.surfaceAreas();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(state.temperature()), "K" });
    table.add_row({ "Pressure", strfix(state.pressure()*1e-5), "bar" });
    table.add_row({ "Charge:", strsci(state.charge()), "mol" });

    table.add_row({ "Element Amount:", "", "" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), strsci(b[i]), "mol" });
    table.add_row({ "Species Amount:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(n[i]), "mol" });

    if(surfaces.size())
    {
        table.add_row({ "Surface Area:", "", "" });
            for(auto [k, pair] : enumerate(surfaces))
                table.add_row({ ":: " + phases[pair.first].name() + " : " + phases[pair.second].name(), strfix(surface_areas[k]), "m2" });
    }

    auto i = 0;
    for(auto& row : table)
    {
        if(i >= 2)  // apply from the third row
            table[i]
                .format()
                .border_top("")
                .column_separator("")
                .corner_top_left("")
                .corner_top_right("");
        i += 1;
    }

    table.row(0).format().font_style({FontStyle::bold});  // Bold face for header
    table.column(1).format().font_align(FontAlign::right); // Value column with right alignment
    table.column(2).format().font_align(FontAlign::right); // Unit column with right alignment

    auto old_locale = std::locale::global(std::locale("C")); // This locale logic is needed to avoid UnicodeDecodeError: 'utf-8' codec can't decode byte 0xa0 in position ...
    out << table;
    std::locale::global(old_locale);

    return out;
}

} // namespace Reaktoro
