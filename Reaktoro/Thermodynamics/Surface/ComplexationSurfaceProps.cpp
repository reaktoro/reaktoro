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
    const auto exchange_phases = system.phases().withAggregateState(AggregateState::Adsorbed);
    warning(exchange_phases.size() > 1,
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

struct ComplexationSurfaceProps::Impl
{
    /// The chemical system to which the complexation surface phase belongs.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the complexation surface phase in the system.
    const Index iphase;

    /// The underlying Phase object for the complexation surface phase in the system.
    const Phase phase;

    /// The complexation surface phase as an complexation surface.
    ComplexationSurface complexation_surface;

    /// The state representing the complexation surface.
    ComplexationSurfaceState complexation_surface_state;

    /// The chemical properties of the complexation surface phase.
    ChemicalPropsPhase props;

    /// The formula matrix of the complexation surface species.
    MatrixXd Aex;

    /// The amounts of the species in the complexation surface phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    ArrayXr nex;

    /// The extra properties and data produced during the evaluation of the complexation surface phase activity model.
    Map<String, Any> extra;

    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexComplexationSurfacePhase(system)),
      phase(system.phase(iphase)),
      props(phase),
      complexation_surface(phase.species())
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == complexation_surface.species().size());

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
            complexation_surface = std::any_cast<ComplexationSurface>(extra["ComplexationSurface"]);
        props.update(T, P, n.segment(ifirst, size), extra);
        update(props);
    }

    /// Update the complexation surface properties with given complexation surface phase properties.
    auto update(const ChemicalPropsPhase& phase_props) -> void
    {
        props = phase_props;
        nex = props.speciesAmounts();

        complexation_surface_state = complexation_surface.state(props.temperature(), props.pressure(), props.speciesMoleFractions());

        // Fetch ionic strength from the aqueous solution in contact with the complexation surface
        if(extra["DiffusiveLayerState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& state = std::any_cast<AqueousMixtureState>(extra["DiffusiveLayerState"]);
            complexation_surface_state.potential(state.Ie);

            // Exit the function without evaluating the next if
            return;
        }
        else if(extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& state = std::any_cast<AqueousMixtureState>(extra["AqueousMixtureState"]);
            complexation_surface_state.potential(state.Ie);
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
        return nex / nex.sum();
    }

    /// Return the fraction of an complexation surface species (in eq).
    auto speciesFraction(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return nex[idx] / nex.sum();
    }

    /// Return the complexation surface state.
    auto complexationSurfaceState() const -> ComplexationSurfaceState
    {
        return complexation_surface_state;
    }

};

ComplexationSurfaceProps::ComplexationSurfaceProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ComplexationSurfaceProps::ComplexationSurfaceProps(const ChemicalState& state)
: pimpl(new Impl(state.system(), state))
{}

ComplexationSurfaceProps::ComplexationSurfaceProps(const ChemicalProps& props)
: pimpl(new Impl(props.system(), props))
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

auto ComplexationSurfaceProps::complexationSurfaceState() const -> ComplexationSurfaceState
    {
        return pimpl->complexationSurfaceState();
    }

auto ComplexationSurfaceProps::phase() const -> const Phase&
{
    return pimpl->phase;
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

    // Auxiliary variables
    const auto F = faradayConstant;
    const auto R = universalGasConstant;

    // Check if the size of the species' and elements' data containers are the same
    assert(species.size() == ns.size());
    assert(elements.size() == ne.size());

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ ":: As    (specific area)" , str(props.complexationSurfaceState().As)       , "m2/kg" });
    table.add_row({ ":: mass  (mass)"          , str(props.complexationSurfaceState().mass)       , "kg" });
    table.add_row({ ":: Z     (charge)"        , str(props.complexationSurfaceState().Z)        , "eq"    });
    table.add_row({ ":: sigma (charge density)", str(props.complexationSurfaceState().sigma)    , "C/m2"  });
    table.add_row({ ":: psi   (potential) "    , str(props.complexationSurfaceState().psi), "Volt"  });
    table.add_row({ ":: -F*psi/RT"             , str(- F * props.complexationSurfaceState().psi / R / props.complexationSurfaceState().T), ""  });
    table.add_row({ ":: exp(-F*psi/RT)"        , str(exp(- F * props.complexationSurfaceState().psi / R / props.complexationSurfaceState().T)), ""  });

    table.add_row({ ":: mass  (mass)" , str(props.complexationSurfaceState().mass)       , "kg" });
    table.add_row({ "Element Amounts:" });
    for(auto i = 0; i < elements.size(); ++i)
        table.add_row({ ":: " + elements[i].symbol(), str(ne[i]), "mole" });
    table.add_row({ "Species Amounts:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].name(), str(ns[i]), "mole" });
    table.add_row({ "Fractions:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].name(), str(x[i]), "" });
    table.add_row({ "Log (base 10) Amount:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].name(), str(log10(ns[i])), "molal" });

    // TODO: figure out how surface complexation species are considered and how is molality calculated?
    //  Site amount:
    // 	  2.500e-05  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           1.147e-05       0.459   1.147e-05      -4.940
    //	Hfo_sOHCa+2       1.005e-05       0.402   1.005e-05      -4.998
    //	Hfo_sOH2+         1.754e-06       0.070   1.754e-06      -5.756
    //	Hfo_sO-           1.719e-06       0.069   1.719e-06      -5.765
    //	Hfo_sOCd+         9.469e-10       0.000   9.469e-10      -9.024

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
