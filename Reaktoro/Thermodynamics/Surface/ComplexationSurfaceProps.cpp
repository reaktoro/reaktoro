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

#include "ComplexationSurfaceProps.hpp"

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
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {
namespace {

/// Return the index of the first complexation surface phase in the system.
auto indexComplexationSurfacePhase(const ChemicalSystem& system) -> Index
{
    const auto sorption_phases = system.phases().withAggregateState(AggregateState::Adsorbed);
    warning(sorption_phases.size() > 1,
            "While creating an ComplexationSurfaceProps object, it has been detected ",
            "more than one complexation surface phase in the system. The ComplexationSurfaceProps object "
            "created will correspond to the first complexation surface phase found.");
    const auto idx = system.phases().findWithAggregateState(AggregateState::Adsorbed);
    error(idx >= system.phases().size(),
          "Could not create an ComplexationSurfaceProps object because there is no "
          "phase in the system with aggregate state value AggregateState::Aqueous.");
    return idx;
}

} // namespace

// Auxiliary constants
const auto F = faradayConstant;
const auto R = universalGasConstant;

struct ComplexationSurfaceProps::Impl
{
    /// The chemical system to which the complexation surface phase belongs.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the complexation surface phase in the system.
    const Index iphase;

    /// The underlying Phase object for the complexation surface phase in the system.
    const Phase phase;

    /// The complexation surface phase as an complexation surface.
    ComplexationSurface surface;

    /// The state representing the complexation surface.
    ComplexationSurfaceState surface_state;

    /// The chemical properties of the complexation surface phase.
    ChemicalPropsPhase props;

    /// The formula matrix of the complexation surface species.
    MatrixXd Aex;

    /// The amounts of the species in the complexation surface phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    ArrayXr nex;

    /// The extra properties and data produced during the evaluation of the complexation surface phase activity model.
    Map<String, Any> extra;

    Impl(const ComplexationSurface& surface, const ChemicalSystem& system)
    : system(system),
      iphase(indexComplexationSurfacePhase(system)),
      phase(system.phase(iphase)),
      props(phase),
      surface(surface)
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == surface.species().size());

        Aex = detail::assembleFormulaMatrix(phase.species(), phase.elements());
    }

    Impl(const ComplexationSurface& surface, const ChemicalSystem& system, const ChemicalState& state)
    : Impl(surface, system)
    {
        update(state);
    }

    Impl(const ComplexationSurface& surface, const ChemicalSystem& system, const ChemicalProps& props)
    : Impl(surface, system)
    {
        update(props);
    }

    /// Update the complexation surface properties with given chemical state.
    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        extra = state.props().extra();

        if(extra["ComplexationSurface"].has_value())
        {
            surface = std::any_cast<ComplexationSurface>(extra["ComplexationSurface"]);
            surface_state = std::any_cast<ComplexationSurfaceState>(extra["ComplexationSurfaceState"]);
        }

        props.update(T, P, n.segment(ifirst, size), extra);
        update(props);
    }

    /// Update the complexation surface properties with given complexation surface phase properties.
    auto update(const ChemicalPropsPhase& phase_props) -> void
    {
        props = phase_props;
        nex = props.speciesAmounts();

        if(extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& state = std::any_cast<AqueousMixtureState>(extra["AqueousMixtureState"]);
            surface_state.updatePotential(state.Ie);
        }
    }

    /// Update the complexation surface properties with given chemical properties of the system.
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

    /// Return the amounts of the species on the complexation surface (in moles).
    auto speciesAmounts() const -> ArrayXr
    {
        return nex;
    }

    /// Return the amounts of an complexation surface species (in moles).
    auto speciesAmount(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return nex[idx];
    }

    /// Return the fraction of the species on the complexation surface composition (in eq).
    auto speciesFractions() const -> ArrayXr
    {
        return surface.sites().size() * (nex / nex.sum());
    }

    /// Return the fraction of an complexation surface species (in eq).
    auto speciesFraction(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return surface.sites().size() * nex[idx] / nex.sum();
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

    /// Return the complexation surface charge.
    auto Z() -> real
    {
        return (surface.charges()*nex).sum();
    }

    /// Return the complexation surface charge density.
    auto charge(real Z) const -> real
    {
        return F*Z/surface.specificSurfaceArea()/surface.mass();
    }

    /// Return the surface complexation potential for given ionic strength of the neighboring phase
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

ComplexationSurfaceProps::ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalSystem& system)
: pimpl(new Impl(surface, system))
{}

ComplexationSurfaceProps::ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalState& state)
: pimpl(new Impl(surface, state.system(), state))
{}

ComplexationSurfaceProps::ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalProps& props)
: pimpl(new Impl(surface, props.system(), props))
{}

ComplexationSurfaceProps::ComplexationSurfaceProps(const ComplexationSurfaceProps& other)
: pimpl(new Impl(*other.pimpl))
{}

ComplexationSurfaceProps::~ComplexationSurfaceProps()
{}

auto ComplexationSurfaceProps::operator=(ComplexationSurfaceProps other) -> ComplexationSurfaceProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ComplexationSurfaceProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto ComplexationSurfaceProps::update(const ChemicalProps& props) -> void
{
    pimpl->update(props);
}

auto ComplexationSurfaceProps::elementAmount(const StringOrIndex& symbol) const -> real
{
    return pimpl->elementAmount(symbol);
}

auto ComplexationSurfaceProps::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto ComplexationSurfaceProps::speciesAmounts() const -> ArrayXr
{
    return pimpl->speciesAmounts();
}

auto ComplexationSurfaceProps::speciesAmount(const StringOrIndex& name) const -> real
{
    return pimpl->speciesAmount(name);
}

auto ComplexationSurfaceProps::speciesFractions() const -> ArrayXr
{
    return pimpl->speciesFractions();
}

auto ComplexationSurfaceProps::speciesFraction(const StringOrIndex& name) const -> real
{
    return pimpl->speciesFraction(name);
}

auto ComplexationSurfaceProps::speciesActivityLg(const StringOrIndex& name) const -> real
{
    return pimpl->speciesActivityLg(name);
}

auto ComplexationSurfaceProps::speciesActivitiesLg() const -> ArrayXr
{
    return pimpl->speciesActivitiesLg();
}

auto ComplexationSurfaceProps::complexationSurfaceState() const -> ComplexationSurfaceState
{
    return pimpl->surface_state;
}

auto ComplexationSurfaceProps::complexationSurface() const -> ComplexationSurface
{
    return pimpl->surface;
}

auto ComplexationSurfaceProps::Z() const -> real
{
    return pimpl->Z();
}

auto ComplexationSurfaceProps::charge(real Z) const -> real
{
    return pimpl->charge(Z);
}

auto ComplexationSurfaceProps::potential(real I, real sigma) const -> real
{
    return pimpl->potential(I, sigma);
}

auto ComplexationSurfaceProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

auto ComplexationSurfaceProps::extra() const -> const Map<String, Any>&
{
    return pimpl->extra;
}

auto ComplexationSurfaceProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ComplexationSurfaceProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const ComplexationSurfaceProps& props) -> std::ostream&
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
    const auto T = props.complexationSurfaceState().T;
    const auto Z = props.Z();
    auto I = 0.0;
    // Export aqueous mixture state via `extra` data member
    if (props.extra().at("AqueousMixtureState").has_value())
        I = std::any_cast<AqueousMixtureState>(props.extra().at("AqueousMixtureState")).Ie;

    // Update charge and potential with Z calculate via amount of species
    const auto sigma = props.charge(Z);
    const auto psi = props.potential(I, sigma);

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ ":: Z     (total charge)"   , str(props.Z()), "eq" });
    table.add_row({ ":: sigma (charge)"        , str(sigma), "C/m2" });
    table.add_row({ ":: psi   (potential) "    , str(psi), "Volt" });
    table.add_row({ ":: ::     -F*psi/(R*T)"   , str(- F*psi/R/T), "" });
    table.add_row({ ":: :: exp(-F*psi/(R*T))"  , str(exp(- F*psi/R/T)), "" });
    table.add_row({ ":: As    (specific area)" , str(props.complexationSurfaceState().As), "m2/kg" });
    table.add_row({ ":: mass  (mass)"          , str(props.complexationSurfaceState().mass), "kg" });

    table.add_row({ "Element Amounts:" });
    for(auto i = 0; i < elements.size(); ++i)
        table.add_row({ ":: " + elements[i].symbol(), str(ne[i]), "mole" });
    table.add_row({ "Species Amounts:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), str(ns[i]), "mole" });
    table.add_row({ "Fractions:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), str(x[i]), "" });

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
