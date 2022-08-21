// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022x Allan Leal
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

#include "ComplexationSurfaceSiteProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {
namespace {

/// Return the index of the first complexation surface site phase in the system.
auto indexComplexationSurfaceSitePhase(const ChemicalSystem& system, const ComplexationSurfaceSite& site) -> Index
{
    const auto exchange_phases = system.phases().withAggregateState(AggregateState::Adsorbed).withNames(site.name());
    error(exchange_phases.size() > 1,
            "While creating an ComplexationSurfaceSiteProps object, it has been detected ",
            "more than one complexation surface site phase with the name " + site.name() +
            " in the system.");
    const auto idx = system.phases().findWithName(site.name());
    error(idx >= system.phases().size(),
          "Could not create an ComplexationSurfaceSiteProps object because there is no "
          "phase in the system with the name " + site.name());
    return idx;
}

} // namespace

struct ComplexationSurfaceSiteProps::Impl
{
    /// The chemical system to which the complexation surface site phase belongs.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the complexation surface site phase in the system.
    const Index iphase;

    /// The underlying Phase object for the complexation surface site phase in the system.
    const Phase phase;

    /// The complexation surface phase as an complexation surface.
    ComplexationSurface surface;

    /// The complexation surface site phase as an complexation surface site.
    ComplexationSurfaceSite site;

    /// The state representing the complexation surface.
    ComplexationSurfaceSiteState site_state;

    /// The chemical properties of the complexation surface site phase.
    ChemicalPropsPhase props;

    /// The formula matrix of the complexation surface site species.
    MatrixXd Aex;

    /// The amounts of the species in the complexation surface site phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    ArrayXr nex;

    /// The extra properties and data produced during the evaluation of the complexation surface site phase activity model.
    Map<String, Any> extra;

    Impl(const ComplexationSurfaceSite& site, const ChemicalSystem& system)
    : system(system),
      iphase(indexComplexationSurfaceSitePhase(system, site)),
      phase(system.phase(iphase)),
      props(phase),
      site(site)
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == site.species().size());

        Aex = detail::assembleFormulaMatrix(phase.species(), phase.elements());
    }

    Impl(const ComplexationSurfaceSite& site, const ChemicalSystem& system, const ChemicalState& state)
    : Impl(site, system)
    {
        update(state);
    }

    Impl(const ComplexationSurfaceSite& site, const ChemicalSystem& system, const ChemicalProps& props)
    : Impl(site, system)
    {
        update(props);
    }

    /// Update the complexation surface site properties with given chemical state.
    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        extra = state.props().extra();

        if(extra["ComplexationSurfaceSiteState" + site.name()].has_value())
            site_state = std::any_cast<ComplexationSurfaceSiteState>(extra["ComplexationSurfaceSiteState" + site.name()]);
        if(extra["ComplexationSurfaceSite" + site.name()].has_value())
            site = std::any_cast<ComplexationSurfaceSite>(extra["ComplexationSurfaceSite" + site.name()]);

        props.update(T, P, n.segment(ifirst, size), extra);
        update(props);
    }

    /// Update the complexation surface site properties with given complexation surface site phase properties.
    auto update(const ChemicalPropsPhase& phase_props) -> void
    {
        props = phase_props;
        nex = props.speciesAmounts();
    }

    /// Update the complexation surface site properties with given chemical properties of the system.
    auto update(const ChemicalProps& sysprops) -> void
    {
        // Copy extra data from the chemical properties class
        extra = sysprops.extra();

        // Update content of the
        update(sysprops.phaseProps(iphase));
    }

    /// Return the amount of an element (in moles).
    auto elementAmount(const StringOrIndex& symbol) const -> real
    {
        const auto idx = detail::resolveElementIndex(phase, symbol);
        return Aex.row(idx) * VectorXr(nex);
    }

    /// Return the amounts of the elements (in moles).
    auto elementAmounts() const -> ArrayXr
    {
        const auto E = phase.elements().size();
        return (Aex.topRows(E) * VectorXr(nex)).array();
    }

    /// Return the amounts of the species on the complexation surface site (in moles).
    auto speciesAmounts() const -> ArrayXr
    {
        return nex;
    }

    /// Return the amounts of an complexation surface site species (in moles).
    auto speciesAmount(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return nex[idx];
    }

    /// Return the fraction of the species on the complexation surface site composition (in eq).
    auto speciesFractions() const -> ArrayXr
    {
        return nex / nex.sum();
    }

    /// Return the fraction of an complexation surface site species (in eq).
    auto speciesFraction(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return nex[idx] / nex.sum();
    }

    /// Return the base-10 logarithm of the activity of the species on the complexation surface.
    auto speciesActivitiesLg() const -> ArrayXr
    {
        auto lng = props.speciesActivitiesLn();
        return lng / std::log(10);
    }

    /// Return the base-10 logarithm of the activity of an surface complexation species.
    auto speciesActivityLg(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        auto lng = props.speciesActivitiesLn()[idx];
        return lng / std::log(10);
    }

    /// Return the complexation surface site charge.
    auto charge() -> real
    {
        return (site.charges()*nex).sum();
    }

    /// Return the complexation surface site charge density.
    auto sigma(real Z) const -> real
    {
        return F*Z/site.specificSurfaceArea()/site.mass();
    }
};

ComplexationSurfaceSiteProps::ComplexationSurfaceSiteProps(const ComplexationSurfaceSite& site, const ChemicalSystem& system)
: pimpl(new Impl(site, system))
{}

ComplexationSurfaceSiteProps::ComplexationSurfaceSiteProps(const ComplexationSurfaceSite& site, const ChemicalState& state)
: pimpl(new Impl(site, state.system(), state))
{}

ComplexationSurfaceSiteProps::ComplexationSurfaceSiteProps(const ComplexationSurfaceSite& site, const ChemicalProps& props)
: pimpl(new Impl(site, props.system(), props))
{}

ComplexationSurfaceSiteProps::ComplexationSurfaceSiteProps(const ComplexationSurfaceSiteProps& other)
: pimpl(new Impl(*other.pimpl))
{}

ComplexationSurfaceSiteProps::~ComplexationSurfaceSiteProps()
{}

auto ComplexationSurfaceSiteProps::operator=(ComplexationSurfaceSiteProps other) -> ComplexationSurfaceSiteProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ComplexationSurfaceSiteProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto ComplexationSurfaceSiteProps::update(const ChemicalProps& props) -> void
{
    pimpl->update(props);
}

auto ComplexationSurfaceSiteProps::elementAmount(const StringOrIndex& symbol) const -> real
{
    return pimpl->elementAmount(symbol);
}

auto ComplexationSurfaceSiteProps::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto ComplexationSurfaceSiteProps::speciesAmounts() const -> ArrayXr
{
    return pimpl->speciesAmounts();
}

auto ComplexationSurfaceSiteProps::speciesAmount(const StringOrIndex& name) const -> real
{
    return pimpl->speciesAmount(name);
}

auto ComplexationSurfaceSiteProps::speciesFractions() const -> ArrayXr
{
    return pimpl->speciesFractions();
}

auto ComplexationSurfaceSiteProps::speciesFraction(const StringOrIndex& name) const -> real
{
    return pimpl->speciesFraction(name);
}

auto ComplexationSurfaceSiteProps::speciesActivityLg(const StringOrIndex& name) const -> real
{
    return pimpl->speciesActivityLg(name);
}

auto ComplexationSurfaceSiteProps::speciesActivitiesLg() const -> ArrayXr
{
    return pimpl->speciesActivitiesLg();
}

auto ComplexationSurfaceSiteProps::complexationSurfaceSiteState() const -> ComplexationSurfaceSiteState
{
    return pimpl->site_state;
}

auto ComplexationSurfaceSiteProps::complexationSurfaceSite() const -> ComplexationSurfaceSite
{
    return pimpl->site;
}

auto ComplexationSurfaceSiteProps::charge() const -> real
{
    return pimpl->charge();
}

auto ComplexationSurfaceSiteProps::sigma(real Z) const -> real
{
    return pimpl->sigma(Z);
}

auto ComplexationSurfaceSiteProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

auto ComplexationSurfaceSiteProps::extra() const -> const Map<String, Any>&
{
    return pimpl->extra;
}

auto ComplexationSurfaceSiteProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ComplexationSurfaceSiteProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const ComplexationSurfaceSiteProps& props) -> std::ostream&
{
    // Extract the output data
    const auto elements = props.phase().elements();
    const auto species = props.phase().species();
    const auto ne = props.elementAmounts();
    const auto ns = props.speciesAmounts();
    const auto x = props.speciesFractions();
    const auto log10a = props.speciesActivitiesLg();

    // Check if the size of the species' and elements' data containers are the same
    assert(species.size() == ns.size());
    assert(elements.size() == ne.size());

    // Fetch auxiliary properties
    const auto T = props.complexationSurfaceSiteState().T;
    const auto charge = props.charge();
    const auto sigma = props.sigma(charge);

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ ":: Z     (total charge)"   , str(charge), "eq" });
    table.add_row({ ":: sigma (charge)"        , str(sigma), "C/m2" });

//    table.add_row({ "Element Amounts:" });
//    for(auto i = 0; i < elements.size(); ++i)
//        table.add_row({ ":: " + elements[i].symbol(), str(ne[i]), "mole" });
    table.add_row({ "Species Amounts:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), str(ns[i]), "mole" });
    table.add_row({ "Fractions:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), str(x[i]), "" });
    table.add_row({ "Log Base 10 of Species Amounts:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), str(log(ns[i])), "" });
    auto i = 0;
    for(auto& row : table)
    {
        if(i >= 2)  // apply from the third row
            table[i].format()
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
