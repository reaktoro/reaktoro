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
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
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

    /// Construct a ChemicalState::Impl instance with given chemical system.
    Impl(const ChemicalSystem& system)
    : system(system), equilibrium(system), props(system)
    {
        n.setConstant(system.species().size(), 1e-16); // set small positive value for initial species amounts
    }

    auto temperature(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive temperature value, ", val, "K, in a ChemicalState object.");
        T = val;
    }

    auto temperature(real val, String unit) -> void
    {
        temperature(units::convert(val, unit, "K"));
    }

    auto pressure(real val) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure value, ", val, "Pa, in a ChemicalState object.");
        P = val;
    }

    auto pressure(real val, String unit) -> void
    {
        error(val <= 0.0, "Cannot set a non-positive pressure "
            "value, ", val, unit, ", in a ChemicalState object.");
        pressure(units::convert(val, unit, "Pa"));
    }

    auto add(String species, real value, String unit) -> void
    {
        const auto ispecies = system.species().index(species);
        add(ispecies, value, unit);
    }

    auto add(Index ispecies, real value, String unit) -> void
    {
        const auto size = system.species().size();
        errorif(ispecies >= size, "Given species index (", ispecies, ") is out-of-bounds (number of species is ", size, ").");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] += amount;
    }

    auto set(String species, real value, String unit) -> void
    {
        const auto ispecies = system.species().index(species);
        set(ispecies, value, unit);
    }

    auto set(Index ispecies, real value, String unit) -> void
    {
        const auto size = system.species().size();
        errorif(ispecies >= size, "Given species index (", ispecies, ") is out-of-bounds (number of species is ", size, ").");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] = amount;
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

    auto setSpeciesAmount(StringOrIndex species, real amount, String unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount);
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species));
        n[ispecies] = units::convert(amount, unit, "mol");
    }

    auto setSpeciesMass(StringOrIndex species, real mass, String unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass);
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species));
        n[ispecies] = units::convert(mass, unit, "kg") / system.species(ispecies).molarMass();
    }

    auto speciesAmountsInPhase(StringOrIndex phase) const -> ArrayXrConstRef
    {
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        return n.segment(start, size);
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

    auto speciesAmount(StringOrIndex species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species));
        return n[ispecies];
    }

    auto speciesMass(StringOrIndex species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndex(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", detail::stringfy(species));
        return n[ispecies] * system.species(ispecies).molarMass();
    }

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

    auto scaleSpeciesAmountsInPhase(StringOrIndex phase, double scalar) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        n.segment(start, size) *= scalar;
    }

    auto scaleVolume(real volume, String unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto current_volume = props.volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseVolume(StringOrIndex phase, real volume, String unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume);
        volume = units::convert(volume, unit, "m3");
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        props.update(T, P, n);
        const auto current_volume = props.phaseProps(iphase).volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidVolume(real volume, String unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_volume =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        const auto& factor = current_fluid_volume > 0.0 ? volume / current_fluid_volume : real(0.0);
        const auto& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        scaleSpeciesAmounts(factor, ifluidspecies);
    }

    auto scaleSolidVolume(real volume, String unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_volume =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        const auto& factor = current_solid_volume > 0.0 ? volume / current_solid_volume : real(0.0);
        const auto& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        scaleSpeciesAmounts(factor, isolidspecies);
    }

    auto scaleMass(real mass, String unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto current_mass = props.mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseMass(StringOrIndex phase, real mass, String unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass);
        mass = units::convert(mass, unit, "kg");
        const auto iphase = detail::resolvePhaseIndex(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", detail::stringfy(phase));
        props.update(T, P, n);
        const auto current_mass = props.phaseProps(iphase).mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidMass(real mass, String unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_mass =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        const auto& factor = current_fluid_mass > 0.0 ? mass / current_fluid_mass : real(0.0);
        const auto& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        scaleSpeciesAmounts(factor, ifluidspecies);
    }

    auto scaleSolidMass(real mass, String unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_mass =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        const auto& factor = current_solid_mass > 0.0 ? mass / current_solid_mass : real(0.0);
        const auto& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        scaleSpeciesAmounts(factor, isolidspecies);
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

auto ChemicalState::temperature(real value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::temperature(real value, String unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::pressure(real value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::pressure(real value, String unit) -> void
{
    pimpl->pressure(value, unit);
}

auto ChemicalState::add(String species, real value, String unit) -> void
{
    pimpl->add(species, value, unit);
}

auto ChemicalState::add(Index ispecies, real value, String unit) -> void
{
    pimpl->add(ispecies, value, unit);
}

auto ChemicalState::set(String species, real value, String unit) -> void
{
    pimpl->set(species, value, unit);
}

auto ChemicalState::set(Index ispecies, real value, String unit) -> void
{
    pimpl->set(ispecies, value, unit);
}

auto ChemicalState::setTemperature(real value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::setTemperature(real value, String unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::setPressure(real value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::setPressure(real value, String unit) -> void
{
    pimpl->pressure(value, unit);
}

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

auto ChemicalState::setSpeciesAmount(StringOrIndex species, real amount, String unit) -> void
{
    pimpl->setSpeciesAmount(species, amount, unit);
}

auto ChemicalState::setSpeciesMass(StringOrIndex species, real mass, String unit) -> void
{
    pimpl->setSpeciesMass(species, mass, unit);
}

auto ChemicalState::scaleSpeciesAmounts(real scalar) -> void
{
    pimpl->scaleSpeciesAmounts(scalar);
}

auto ChemicalState::scaleSpeciesAmountsInPhase(StringOrIndex phase, real scalar) -> void
{
    pimpl->scaleSpeciesAmountsInPhase(phase, scalar);
}

auto ChemicalState::scaleVolume(real value, String unit) -> void
{
    pimpl->scaleVolume(value, unit);
}

auto ChemicalState::scalePhaseVolume(StringOrIndex phase, real value, String unit) -> void
{
    pimpl->scalePhaseVolume(phase, value, unit);
}

auto ChemicalState::scaleFluidVolume(real value, String unit) -> void
{
    pimpl->scaleFluidVolume(value, unit);
}

auto ChemicalState::scaleSolidVolume(real value, String unit) -> void
{
    pimpl->scaleSolidVolume(value, unit);
}

auto ChemicalState::scaleMass(real value, String unit) -> void
{
    pimpl->scaleMass(value, unit);
}

auto ChemicalState::scalePhaseMass(StringOrIndex phase, real value, String unit) -> void
{
    pimpl->scalePhaseMass(phase, value, unit);
}

auto ChemicalState::scaleFluidMass(real value, String unit) -> void
{
    pimpl->scaleFluidMass(value, unit);
}

auto ChemicalState::scaleSolidMass(real value, String unit) -> void
{
    pimpl->scaleSolidMass(value, unit);
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

auto ChemicalState::speciesAmountsInPhase(StringOrIndex phase) const -> ArrayXrConstRef
{
    return pimpl->speciesAmountsInPhase(phase);
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

auto ChemicalState::speciesAmount(StringOrIndex species) const -> real
{
    return pimpl->speciesAmount(species);
}

auto ChemicalState::speciesMass(StringOrIndex species) const -> real
{
    return pimpl->speciesMass(species);
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
    const auto n = state.speciesAmounts();
    const auto b = state.elementAmounts();
    const auto species = state.system().species();
    const auto elements = state.system().elements();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", str(state.temperature()), "K" });
    table.add_row({ "Pressure", str(state.pressure()), "Pa" });
    table.add_row({ "Charge:", str(state.charge()), "mol" });

    table.add_row({ "Element Amount:", "", "" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), str(b[i]), "mol" });
    table.add_row({ "Species Amount:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(n[i]), "mol" });

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
