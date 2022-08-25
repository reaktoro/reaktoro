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

#include "SurfaceProps.hpp"

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
#include <Reaktoro/Thermodynamics/Surface/Surface.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {

namespace {

/// Return the index of the first surface site phase in the system.
auto indexSurfaceSitePhase(const ChemicalSystem& system, const SurfaceSite& site) -> Index
{
    const auto exchange_phases = system.phases().withAggregateState(AggregateState::Adsorbed).withNames(site.name());
    error(exchange_phases.size() > 1,
          "While creating an SurfaceProps object, it has been detected that",
          " more than one surface site phase with the name " + site.name() +
          " in the system.");
    const auto idx = system.phases().findWithName(site.name());
    error(idx >= system.phases().size(),
          "Could not create an SurfaceProps object because there is no "
          "phase in the system with the name " + site.name());
    return idx;
}

// Function formatting table's header and columns' alignment
auto formatTable(Table& table) -> void
{
    auto i = 0;
    for (auto &row : table) {
        if (i >= 2)  // apply from the third row
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
}

} //

struct SurfaceProps::Impl
{
    /// The chemical system to which the surface phase belongs.
    const ChemicalSystem system;

    /// The indices of the underlying Phases object for the surface site phases in the system.
    std::map<std::string, Index> iphases;

    /// The underlying Phase object for the surface phase in the system.
    std::map<std::string, Phase> phases;

    /// The surface.
    Surface surface;

    /// The state representing the surface.
    SurfaceState surface_state;

    /// The chemical properties of the system.
    ChemicalProps props;

    /// The formula matrix of the surface species.
    std::map<std::string, MatrixXd> Aex;

    /// The amounts of the species in the surface phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    std::map<std::string, ArrayXr> nex;

    /// The extra properties and data produced during the evaluation of the surface phase activity model.
    Map<String, Any> extra;

    Impl(const Surface& surface, const ChemicalSystem& system)
    : system(system),
      props(system),
      surface(surface)
    {

        for(auto [tag, site] : surface.sites())
        {
            auto iphase_site = indexSurfaceSitePhase(system, site);
            auto phase_site = Phase(system.phase(iphase_site));
            auto props_site = props.phaseProps(iphase_site);

            assert(phase_site.species().size() > 0);
            assert(phase_site.species().size() == props_site.phase().species().size());
            assert(phase_site.species().size() == site.species().size());

            Aex[tag] = detail::assembleFormulaMatrix(phase_site.species(), phase_site.elements());

            iphases[tag] = iphase_site;
            phases[tag] = phase_site;
        }
    }

    Impl(const Surface& surface, const ChemicalSystem& system, const ChemicalState& state)
    : Impl(surface, system)
    {
        update(state);
    }

    /// Update the surface properties with given chemical state.
    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();

        for(auto [tag, site] : surface.sites())
        {
            const auto ifirst = system.phases().numSpeciesUntilPhase(iphases[tag]);
            const auto size = phases[tag].species().size();

            auto props_site = props.phaseProps(iphases[tag]);
            props_site.update(T, P, n.segment(ifirst, size), extra);

            nex[tag] = props_site.speciesAmounts();
        }

        extra = state.props().extra();

        if(extra["Surface"].has_value())
        {
            surface = std::any_cast<Surface>(extra["Surface"]);
            surface_state = std::any_cast<SurfaceState>(extra["SurfaceState"]);
        }
    }

    /// Return the amount of an element (in moles).
    auto elementAmount(const String& site_tag, const StringOrIndex& symbol) const -> real
    {
        const auto idx = detail::resolveElementIndex(phases.at(site_tag), symbol);
        return Aex.at(site_tag).row(idx) * VectorXr(nex.at(site_tag));
    }

    /// Return the amounts of the elements (in moles).
    auto elementAmounts(const String& site_tag) const -> ArrayXr
    {
        const auto E = phases.at(site_tag).elements().size();
        return (Aex.at(site_tag).topRows(E) * VectorXr(nex.at(site_tag))).array();
    }

    /// Return the amounts of the species on the surface (in moles).
    auto speciesAmounts(const String& site_tag) const -> ArrayXr
    {
        return nex.at(site_tag);
    }

    /// Return the amounts of an surface species (in moles).
    auto speciesAmount(const String& site_tag, const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phases.at(site_tag), name);
        return nex.at(site_tag)[idx];
    }

    /// Return the fraction of the species on the surface composition (in eq).
    auto speciesFractions(const String& site_tag) const -> ArrayXr
    {
        return nex.at(site_tag) / nex.at(site_tag).sum();
    }

    /// Return the fraction of an surface species (in eq).
    auto speciesFraction(const String& site_tag, const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phases.at(site_tag), name);
        return nex.at(site_tag)[idx] / nex.at(site_tag).sum();
    }

    /// Return the surface charge.
    auto Z() -> real
    {
        real Z = 0;
        for(auto [tag, site] : surface.sites())
            Z += (site.charges()*nex.at(tag)).sum();
        return Z;
    }

    /// Return the surface charge density.
    auto charge(real Z) const -> real
    {
        return F*Z/surface.specificSurfaceArea()/surface.mass();
    }

    /// Return the surface potential for given ionic strength of the neighboring phase
    auto potential(real I, real sigma) const -> real
    {
        // Auxiliary variables
        const auto T = surface_state.T;

        // Using formula sigma = 0.1174*I^0.5*sinh(F*potential/R/T/2) and arcsinh(y) = ln(y+(y^2+1)^1⁄2)
        const auto y = sigma/(0.1174*sqrt(I));
        const auto arcsinhy = asinh(y);
        return 2*R*T*arcsinhy/F;
    }
};

SurfaceProps::SurfaceProps(const Surface& surface, const ChemicalSystem& system)
: pimpl(new Impl(surface, system))
{}

SurfaceProps::SurfaceProps(const Surface& surface, const ChemicalState& state)
: pimpl(new Impl(surface, state.system(), state))
{}

SurfaceProps::SurfaceProps(const SurfaceProps& other)
: pimpl(new Impl(*other.pimpl))
{}

SurfaceProps::~SurfaceProps()
{}

auto SurfaceProps::operator=(SurfaceProps other) -> SurfaceProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SurfaceProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto SurfaceProps::elementAmount(const String& site_tag, const StringOrIndex& symbol) const -> real
{
    return pimpl->elementAmount(site_tag, symbol);
}

auto SurfaceProps::elementAmounts(const String& site_tag) const -> ArrayXr
{
    return pimpl->elementAmounts(site_tag);
}

auto SurfaceProps::speciesAmounts(const String& site_tag) const -> ArrayXr
{
    return pimpl->speciesAmounts(site_tag);
}

auto SurfaceProps::speciesAmount(const String& site_tag, const StringOrIndex& name) const -> real
{
    return pimpl->speciesAmount(site_tag, name);
}

auto SurfaceProps::speciesFractions(const String& site_tag) const -> ArrayXr
{
    return pimpl->speciesFractions(site_tag);
}

auto SurfaceProps::speciesFraction(const String& site_tag, const StringOrIndex& name) const -> real
{
    return pimpl->speciesFraction(site_tag, name);
}

auto SurfaceProps::surfaceState() const -> SurfaceState
{
    return pimpl->surface_state;
}

auto SurfaceProps::surface() const -> Surface
{
    return pimpl->surface;
}

auto SurfaceProps::Z() const -> real
{
    return pimpl->Z();
}

auto SurfaceProps::charge(real Z) const -> real
{
    return pimpl->charge(Z);
}

auto SurfaceProps::potential(real I, real sigma) const -> real
{
    return pimpl->potential(I, sigma);
}

auto SurfaceProps::phase(const String& site_tag) const -> const Phase&
{
    return pimpl->phases.at(site_tag);
}

auto SurfaceProps::extra() const -> const Map<String, Any>&
{
    return pimpl->extra;
}

auto SurfaceProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto SurfaceProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const SurfaceProps& props) -> std::ostream&
{
    // Fetch auxiliary properties
    auto surface = props.surface();
    auto surface_state = props.surfaceState();

    const auto T = surface_state.T;
    const auto Z = props.Z();
    auto I = 0.0;
    // Export aqueous mixture state via `extra` data member
    if (props.extra().at("AqueousMixtureState").has_value())
        I = std::any_cast<AqueousMixtureState>(props.extra().at("AqueousMixtureState")).Ie;

    // Update charge and potential with Z calculate via amount of species
    const auto sigma = props.charge(Z);
    const auto psi = props.potential(I, sigma);

    Table table_surface;
    table_surface.add_row({ "Properties of the surface " + surface.name(), "Value", "Unit" });
    table_surface.add_row({ ":: Z     (total charge)"   , str(props.Z()), "eq" });
    table_surface.add_row({ ":: sigma (charge)"        , str(sigma), "C/m2" });
    table_surface.add_row({ ":: psi   (potential) "    , str(psi), "Volt" });
    table_surface.add_row({ ":: ::     -F*psi/(R*T)"   , str(- F*psi/R/T), "" });
    table_surface.add_row({ ":: :: exp(-F*psi/(R*T))"  , str(exp(- F*psi/R/T)), "" });
    table_surface.add_row({":: As    (specific area)" , str(props.surfaceState().As), "m2/kg" });
    table_surface.add_row({":: mass  (mass)"          , str(props.surfaceState().mass), "kg" });

    formatTable(table_surface);

    out << table_surface << "\n";

    for(auto [tag, site] : surface.sites()) {

        Table table_site;
        table_site.add_row({ "Properties of the site " + site.name(), "Value", "Unit" });

        // Extract the output data
        const auto elements = props.phase(tag).elements();
        const auto species = props.phase(tag).species();
        const auto ne = props.elementAmounts(tag);
        const auto ns = props.speciesAmounts(tag);
        const auto x = props.speciesFractions(tag);

        // Check if the size of the species' and elements' data containers are the same
        assert(species.size() == ns.size());
        assert(elements.size() == ne.size());

        table_site.add_row({"Element Amounts:"});
        for (auto i = 0; i < elements.size(); ++i)
            table_site.add_row({":: " + elements[i].symbol(), str(ne[i]), "mole"});
        table_site.add_row({"Species Amounts:"});
        for (auto i = 0; i < species.size(); ++i)
            table_site.add_row({":: " + species[i].repr(), str(ns[i]), "mole"});
        table_site.add_row({"Fractions:"});
        for (auto i = 0; i < species.size(); ++i)
            table_site.add_row({":: " + species[i].repr(), str(x[i]), ""});
        table_site.add_row({ "Log Base 10 of Species Amounts:" });
        for(auto i = 0; i < species.size(); ++i)
            table_site.add_row({ ":: " + species[i].repr(), str(log(ns[i])), "" });

        formatTable(table_site);

        out << table_site << "\n";
    }

    return out;
}

} // namespace Reaktoro
