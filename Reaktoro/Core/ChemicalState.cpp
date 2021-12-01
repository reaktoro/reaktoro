// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
        n.setZero(system.species().size());
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

    auto scaleSpeciesAmounts(double scalar) -> void
    {
        Assert(scalar >= 0.0, "Cannot scale the molar amounts of the species.",
               "The given scalar is negative.");
        for(int i = 0; i < n.rows(); ++i)
            setSpeciesAmount(i, speciesAmount(i) * scalar);
    }

    auto scaleSpeciesAmountsInPhase(Index index, double scalar) -> void
    {
        Assert(scalar >= 0.0,
               "Cannot scale the molar amounts of the species.",
               "The given scalar `" + std::to_string(scalar) << "` is negative.")
        Assert(index < system.phases().size(), "Cannot set the volume of the phase.",
               "The given phase index is out of range.")
        const Index start = detail::resolveSpeciesIndex(system, system.phase(index).species()[0].name());
        const Index size = system.phase(index).species().size();
        for(unsigned i = 0; i < size; ++i)
            setSpeciesAmount(start + i, speciesAmount(start + i) * scalar);
    }

    auto scaleVolume(double volume) -> void
    {
        Assert(volume >= 0.0,
               "Cannot set the volume of the chemical state.",
               "The given volume is negative.")
        const auto vtotal = props.volume();
        const auto scalar = (vtotal != 0.0) ? volume/vtotal : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scaleVolume(double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        return scaleVolume(volume);
    }

    auto scalePhaseVolume(Index index, double volume) -> void
    {
        Assert(volume >= 0.0, "Cannot set the volume of the phase.",
               "The given volume is negative.");
        Assert(index < system.phases().size(), "Cannot set the volume of the phase.",
               "The given phase index is out of range.");
        const auto v = props.phaseProps(index).volume();
        const auto scalar = (v != 0.0) ? real(volume)/v : real(0.0);
        scaleSpeciesAmountsInPhase(index, scalar);
    }

    auto scalePhaseVolume(Index index, double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scalePhaseVolume(index, volume);
    }

    auto scalePhaseVolume(std::string name, double volume) -> void
    {
        const Index index = system.indexPhaseWithError(name);
        scalePhaseVolume(index, volume);
    }

    auto scalePhaseVolume(std::string name, double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scalePhaseVolume(name, volume);
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
        const auto& A = system.formulaMatrixElements();
        return A * n.matrix();
    }

    auto elementAmountsInSpecies(const Indices& indices) const -> ArrayXr
    {
        return system.elementAmountsInSpecies(indices, n);
    }

    auto charge() const -> real
    {
        const auto& Az = system.formulaMatrixCharge();
        return (Az * n.matrix())[0];
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

auto ChemicalState::scaleSpeciesAmounts(double scalar) -> void
{
    pimpl->scaleSpeciesAmounts(scalar);
}

auto ChemicalState::scaleSpeciesAmountsInPhase(Index index, double scalar) -> void
{
    pimpl->scaleSpeciesAmountsInPhase(index, scalar);
}

auto ChemicalState::scaleVolume(double volume) -> void
{
    pimpl->scaleVolume(volume);
}

auto ChemicalState::scaleVolume(double volume, std::string units) -> void
{
    pimpl->scaleVolume(volume, units);
}

auto ChemicalState::scalePhaseVolume(Index index, double volume) -> void
{
    pimpl->scalePhaseVolume(index, volume);
}

auto ChemicalState::scalePhaseVolume(Index index, double volume, std::string units) -> void
{
    pimpl->scalePhaseVolume(index, volume, units);
}

auto ChemicalState::scalePhaseVolume(std::string name, double volume) -> void
{
    pimpl->scalePhaseVolume(name, volume);
}

auto ChemicalState::scalePhaseVolume(std::string name, double volume, std::string units) -> void
{
    pimpl->scalePhaseVolume(name, volume, units);
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

auto ChemicalState::elementAmountsInSpecies(const Indices& ispecies) const -> ArrayXr
{
    return pimpl->elementAmountsInSpecies(ispecies);
}

auto ChemicalState::charge() const -> real
{
    return pimpl->charge();
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

auto ChemicalState::props() const -> const ChemicalProps&
{
    return pimpl->props;
}

auto ChemicalState::props() -> ChemicalProps&
{
    return pimpl->props;
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
    pimpl->optstate.p = p;
}

auto ChemicalState::Equilibrium::setControlVariablesQ(ArrayXdConstRef q) -> void
{
    pimpl->q = q;
}

auto ChemicalState::Equilibrium::setOptimaState(const Optima::State& state) -> void
{
    pimpl->optstate = state;
    const auto Nq = state.x.size() - pimpl->Nn;
    pimpl->q = state.x.tail(Nq);
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
    return pimpl->optstate.p;
}

auto ChemicalState::Equilibrium::controlVariablesQ() const -> ArrayXdConstRef
{
    return pimpl->q;
}

auto ChemicalState::Equilibrium::p() const -> ArrayXdConstRef
{
    return pimpl->optstate.p;
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
    const auto T = state.temperature();
    const auto P = state.pressure();
    const auto n = state.speciesAmounts();
    const auto b = state.elementAmounts();
    const auto species = state.system().species();
    const auto elements = state.system().elements();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", str(T), "K" });
    table.add_row({ "Pressure", str(P), "Pa" });
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

    out << table;
    return out;
}

} // namespace Reaktoro
