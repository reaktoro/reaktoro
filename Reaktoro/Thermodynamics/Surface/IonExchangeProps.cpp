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

#include "IonExchangeProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>

namespace Reaktoro {
namespace {

/// Return the index of the first ion exchange phase in the system.
auto indexIonExchangePhase(const ChemicalSystem& system) -> Index
{
    const auto exchange_phases = system.phases().withAggregateState(AggregateState::IonExchange);
    warning(exchange_phases.size() > 1,
        "While creating an IonExchangeProps object, it has been detected ",
        "more than one aqueous phase in the system. The IonExchangeProps object "
        "created will correspond to the first aqueous phase found.");
    const auto idx = system.phases().findWithAggregateState(AggregateState::IonExchange);
    error(idx >= system.phases().size(),
        "Could not create an IonExchangeProps object because there is no "
        "phase in the system with aggregate state value AggregateState::Aqueous.");
    return idx;
}

} // namespace

struct IonExchangeProps::Impl
{
    /// The chemical system to which the ion exchange phase belongs.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the ion exchange phase in the system.
    const Index iphase;

    /// The underlying Phase object for the aqueous phase in the system.
    const Phase phase;

    /// The ion exchange phase as an ion exchange surface.
    IonExchangeSurface exsurface;

    /// The state representing the ion exchange surface.
    IonExchangeSurfaceState exstate;

    /// The chemical properties of the ion exchange phase.
    ChemicalPropsPhase props;

    /// The formula matrix of the aqueous species.
    MatrixXd Aex;

    /// The amounts of the species in the aqueous phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    ArrayXr nex;

    /// The extra properties and data produced during the evaluation of the ion exchange phase activity model.
    Map<String, Any> extra;

    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexIonExchangePhase(system)),
      phase(system.phase(iphase)),
      props(phase),
      exsurface(phase.species())
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == exsurface.species().size());

        Aex = detail::assembleFormulaMatrix(phase.species(), phase.elements());
    }

    Impl(const ChemicalSystem& system, const ChemicalState& state)
    : Impl(system)
    {
        update(state);
    }

    Impl(const ChemicalSystem& system, const ChemicalProps& props)
    : Impl(system)
    {
        update(props);
    }

    /// Update the ion exchange properties with given chemical state.
    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        extra = state.props().extra();
        props.update(T, P, n.segment(ifirst, size), extra);
        update(props);
    }

    /// Update the ion exchange properties with given ion exchange phase properties.
    auto update(const ChemicalPropsPhase& exprops) -> void
    {
        props = exprops;
        exsurface.state(props.temperature(), props.pressure(), props.speciesMoleFractions());
        exstate = exsurface.state();
        nex = props.speciesAmounts();
    }

    /// Update the ion exchange properties with given chemical properties of the system.
    auto update(const ChemicalProps& sysprops) -> void
    {
        // Copy extra data from the chemical properties class
        extra = sysprops.extra();

        // Update content of the
        update(sysprops.phaseProps(iphase));
    }

    /// Return the amount of an element (in moles).
    auto elementAmount(const String& symbol) const -> real
    {
        const auto idx = phase.elements().indexWithSymbol(symbol);
        return Aex.row(idx) * VectorXr(nex);
    }

    /// Return the amounts of the elements (in moles).
    auto elementAmounts() const -> ArrayXr
    {
        const auto E = phase.elements().size();
        return (Aex.topRows(E) * VectorXr(nex)).array();
    }

    /// Return the amounts of the species on the ion exchange surface (in moles).
    auto speciesAmounts() const -> ArrayXr
    {
        return nex;
    }

    /// Return the amounts of an ion exchange species (in moles).
    auto speciesAmount(const String& name) const -> real
    {
        const auto idx = phase.species().indexWithName(name);
        return nex[idx];
    }

    /// Return the equivalences of the species on the ion exchange composition (in eq).
    auto speciesEquivalences() const -> ArrayXr
    {
        // Note: this definition eq = n * ze is consistent with the PHREEQC output,
        // but meq is usually defined via molalities as meq = 1e-3 * m, where m is the molality
        return nex * exsurface.ze();
    }

    /// Return the equivalence of an ion exchange species (in eq).
    auto speciesEquivalence(const String& name) const -> real
    {
        // Note: this definition eq = n * ze is consistent with the PHREEQC output,
        // but meq is usually defined via molalities as meq = 1e-3 * m, where m is the molality
        const auto idx = phase.species().indexWithName(name);
        return nex[idx] * exsurface.ze()[idx];
    }

    /// Return the equivalent fractions of the species on the ion exchange surface.
    auto speciesEquivalentFractions() const -> ArrayXr
    {
        return exstate.beta;
    }

    /// Return the equivalent fraction of an ion exchange species.
    auto speciesEquivalentFraction(const String& name) const -> real
    {
        const auto idx = phase.species().indexWithName(name);
        return exstate.beta[idx];
    }

    /// Return the base-10 logarithm of the activity coefficients of the species on the ion exchange surface.
    auto speciesLog10Gammas() const -> ArrayXr
    {
        auto extr = extra;
        if(!extr.empty() && extr["IonExchangeLnGamma"].has_value())
        {
            // Fetch ln(gamma) from `extra` data-field and convert it to the log10(gamma)
            auto lng = std::any_cast<ArrayXrRef>(extr["IonExchangeLnGamma"]);
            return lng / std::log(10);
        }
        warningif(!extr["IonExchangeLnGamma"].has_value(),
                  "There is no information on the activity coefficients for the ion exchange species, "
                  "so log10(gamma) are set to zero.")
        return ArrayXr::Zero(phase.elements().size());
    }

    /// Return the base-10 logarithm of the activity coefficients of an ion exchange species.
    auto speciesLog10Gamma(const String& name) const -> real
    {
        const auto idx = phase.species().indexWithName(name);
        auto extr = extra;

        if(!extr.empty() && extr["IonExchangeLnGamma"].has_value())
        {
            // Fetch ln(gamma) from `extra` data-field and convert it to the log10(gamma) with the index `idx`
            auto lng = std::any_cast<ArrayXrRef>(extr["IonExchangeLnGamma"])[idx];
            return lng / std::log(10);
        }
        warningif(!extr["IonExchangeLnGamma"].has_value(),
                  "There is no information on the activity coefficients for an ion exchange species " + name +
                  ", so log10(gamma(" + name + " )) is set to zero.")
        return 0;
    }
};

IonExchangeProps::IonExchangeProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

IonExchangeProps::IonExchangeProps(const ChemicalState& state)
: pimpl(new Impl(state.system(), state))
{}

IonExchangeProps::IonExchangeProps(const ChemicalProps& props)
: pimpl(new Impl(props.system(), props))
{}

IonExchangeProps::IonExchangeProps(const IonExchangeProps& other)
: pimpl(new Impl(*other.pimpl))
{}

IonExchangeProps::~IonExchangeProps()
{}

auto IonExchangeProps::operator=(IonExchangeProps other) -> IonExchangeProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto IonExchangeProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto IonExchangeProps::update(const ChemicalProps& props) -> void
{
    pimpl->update(props);
}

auto IonExchangeProps::elementAmount(const String& symbol) const -> real
{
    return pimpl->elementAmount(symbol);
}

auto IonExchangeProps::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto IonExchangeProps::speciesAmounts() const -> ArrayXr
{
    return pimpl->speciesAmounts();
}

auto IonExchangeProps::speciesAmount(const String& name) const -> real
{
    return pimpl->speciesAmount(name);
}

auto IonExchangeProps::speciesEquivalences() const -> ArrayXr
{
    return pimpl->speciesEquivalences();
}

auto IonExchangeProps::speciesEquivalence(const String& name) const -> real
{
    return pimpl->speciesEquivalence(name);
}

auto IonExchangeProps::speciesEquivalentFractions() const -> ArrayXr
{
    return pimpl->speciesEquivalentFractions();
}

auto IonExchangeProps::speciesEquivalentFraction(const String& name) const -> real
{
    return pimpl->speciesEquivalentFraction(name);
}

auto IonExchangeProps::speciesLog10Gamma(const String& name) const -> real
{
    return pimpl->speciesLog10Gamma(name);
}

auto IonExchangeProps::speciesLog10Gammas() const -> ArrayXr
{
    return pimpl->speciesLog10Gammas();
}

auto IonExchangeProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

auto IonExchangeProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto IonExchangeProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const IonExchangeProps& props) -> std::ostream&
{
    // Extract the output data
    const auto elements = props.phase().elements();
    const auto species = props.phase().species();
    const auto ne = props.elementAmounts();
    const auto ns = props.speciesAmounts();
    const auto eq = props.speciesEquivalences();
    const auto beta = props.speciesEquivalentFractions();
    const auto log10g = props.speciesLog10Gammas();

    // Check if the size of the species' and elements' data containers are the same
    assert(species.size() == ns.size());
    assert(elements.size() == ne.size());

    // Create first output table for the elements taking part in the ion exchange
    Table table;
    table.add_row({ "Element", "Amounts [moles]"});
    for(auto i = 0; i < elements.size(); ++i)
        table.add_row({ ":: " + elements[i].symbol(), str(ne[i]), "moles" });

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

    out << table << "\n";

    // Create second output table for the species taking part in the ion exchange
    Table table_species;
    table_species.add_row({"Species" , "Amounts [moles]", "Equivalence [eq]", "Equivalent Fraction", "Log10 Gamma"});
    for(auto i = 0; i < species.size(); ++i)
        table_species.add_row({":: " + species[i].name(), str(ns[i]), str(eq[i]), str(beta[i]), str(log10g[i])});

    auto j = 0;
    for(auto& row : table_species)
    {
        if(j >= 2)  // apply from the third row
            table_species[j]
            .format()
            .border_top("")
            .column_separator("")
            .corner_top_left("")
            .corner_top_right("");
        j += 1;
    }

    table_species.row(0).format().font_style({FontStyle::bold});  // Bold face for header
    table_species.column(1).format().font_align(FontAlign::right); // Amount column with right alignment
    table_species.column(2).format().font_align(FontAlign::right); // Equivalence column with right alignment
    table_species.column(3).format().font_align(FontAlign::right); // Equivalent fraction column with right alignment
    table_species.column(4).format().font_align(FontAlign::right); // Log10 gamma column with right alignment

    out << table_species;

    return out;
}

} // namespace Reaktoro
