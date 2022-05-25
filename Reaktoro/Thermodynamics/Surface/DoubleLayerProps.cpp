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

#include "DoubleLayerProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Optima includes
#include <Optima/Echelonizer.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {
namespace {

/// Return the index of the first double layer phase in the system.
auto indexDoubleLayerPhase(const ChemicalSystem& system) -> Index
{
    const auto idx = system.phases().findWithName("DoubleLayerPhase");
    error(idx >= system.phases().size(),
        "Could not create an DoubleLayerProps object because there is no "
        "phase in the system with a name DoubleLayerPhase.");
    return idx;
}

} // namespace

struct DoubleLayerProps::Impl
{
    /// The chemical system in which the double layer phase is.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the double layer phase in the system.
    const Index iphase;

    /// The underlying Phase object for the double layer phase in the system.
    const Phase phase;

    /// The phase as an aqueous solution.
    const AqueousMixture aqsolution;

    /// The index of the aqueous solvent species H2O in the double layer phase (not in the system!)
    const Index iH2O;

    /// The index of the aqueous solute species H+ in the double layer phase (not in the system!)
    const Index iH;

    /// The chemical properties of the system.
    ChemicalProps props;

    /// The state of the aqueous solution.
    AqueousMixtureState aqstate;

    /// The amounts of the species in the double layer phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    VectorXd naq;

    /// The chemical potentials of the elements in the aqueous phase
    VectorXr lambda;

    /// The non-aqueous species in the database for which saturation indices are calculated.
    SpeciesList nonaqueous;

    /// The formula matrix of the aqueous species in the aqueous phase.
    MatrixXd Aaqs;

    /// The formula matrix of the non-aqueous species for the computation of saturation indices.
    MatrixXd Anon;

    /// The echelon form of the formula matrix `Aaqs` of the aqueous species.
    Optima::Echelonizer echelonizer;

    // The chemical potential models for the non-aqueous species (as if they were pure phases) for the computation of their saturation indices.
    Vec<Fn<real(ChemicalProps const&)>> chemical_potential_models;

    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexDoubleLayerPhase(system)),
      phase(system.phase(iphase)),
      props(system),
      aqsolution(phase.species()),
      iH2O(phase.species().findWithFormula("H2O")),
      iH(phase.species().findWithFormula("H+"))
    {
        const auto size = phase.species().size();

        error(iH2O >= size, "Cannot create DoubleLayerProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H2O.");

        error(iH >= size, "Cannot create DoubleLayerProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H+.");

        // The symbols of the elements in the aqueous phase
        const auto symbols = vectorize(phase.elements(), RKT_LAMBDA(x, x.symbol()));

        // Collect the species from the database that contains the elements in the aqueous phase
        const auto species_same_elements = system.database().species().withElements(symbols);

        // Collect the non-aqueous species from the database that contains the elements in the aqueous phase
        nonaqueous = remove(species_same_elements, RKT_LAMBDA(x, x.aggregateState() == AggregateState::Aqueous));

        // Ensure non-aqueous species are sorted by aggregate state (gases, solids, etc)
        std::sort(nonaqueous.begin(), nonaqueous.end(),
            [](auto l, auto r)
                { return l.aggregateState() < r.aggregateState(); });

        // Assemble the formula matrices of the aqueous and non-aqueous species w.r.t. elements in the aqueous phase
        Aaqs = detail::assembleFormulaMatrix(phase.species(), phase.elements());
        Anon = detail::assembleFormulaMatrix(nonaqueous, phase.elements());

        // Initialize the aqueous state properties
        aqstate.T = NaN;
        aqstate.P = NaN;
        aqstate.rho = NaN;
        aqstate.epsilon = NaN;
        aqstate.Ie = NaN;
        aqstate.Is = NaN;
        aqstate.m.setConstant(size, NaN);
        aqstate.ms.setConstant(size, NaN);
        aqstate.z.setConstant(size, NaN);

        // Compute the initial echelon form of formula matrix `Aaqs`
        echelonizer.compute(Aaqs);
    }

    Impl(ChemicalState const& state)
    : Impl(state.system())
    {
        update(state);
    }

    Impl(ChemicalProps const& props)
    : Impl(props.system())
    {
        update(props);
    }

    auto update(const ChemicalState& state) -> void
    {
        props.update(state);
        update(props);
    }

    auto update(const ChemicalProps& cprops) -> void
    {
        // Auxiliary variables
        auto const& aqprops = cprops.phaseProps(iphase);
        auto const& T = aqprops.temperature();
        auto const& P = aqprops.pressure();
        auto const& x = aqprops.speciesMoleFractions();

        // Update the internal properties of the aqueous phase
        props = cprops;

        // Update the internal aqueous state object
        aqstate = aqsolution.state(T, P, x);

        // Update auxiliary vector naq to be used in the echelonization below
        naq = aqprops.speciesAmounts();

        // Update the echelon form and also the list of basic species
        echelonizer.updateWithPriorityWeights(naq);

        // Compute chemical potentials of the elements in the aqueous phase
        const auto u = aqprops.speciesChemicalPotentials();
        const auto ib = echelonizer.indicesBasicVariables();
        const auto R = echelonizer.R();
        const auto Rb = R.topRows(ib.size());
        const VectorXr ub = u(ib);
        lambda = Rb.transpose() * ub;
    }

    auto temperature() const -> real
    {
        return props.temperature();
    }

    auto pressure() const -> real
    {
        return props.pressure();
    }

    auto elementMolality(const StringOrIndex& symbol) const -> real
    {
        const auto idx = detail::resolveElementIndex(phase, symbol);
        const auto& m = aqstate.m.matrix();
        return Aaqs.row(idx) * m;
    }

    auto elementMolalities() const -> ArrayXr
    {
        const auto E = phase.elements().size();
        const auto& m = aqstate.m.matrix();
        return Aaqs.topRows(E) * m;
    }

    auto speciesMolality(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return aqstate.m[idx];
    }

    auto speciesMolalities() const -> ArrayXr
    {
        return aqstate.m;
    }

    auto speciesCharges() const -> ArrayXr
    {
        return aqstate.z;
    }

    auto ionicStrength() const -> real
    {
        return aqstate.Ie;
    }

    auto ionicStrengthStoichiometric() const -> real
    {
        return aqstate.Is;
    }

    auto pH() const -> real
    {
        auto const& aqprops = props.phaseProps(iphase);
        auto const& ln_aH = aqprops.speciesActivitiesLn()[iH];
        return -ln_aH/ln10;
    }

    auto pE() const -> real
    {
        const auto T = props.temperature();
        const auto E = phase.elements().size();
        const auto lambdaZ = lambda[E];
        const auto RT = universalGasConstant * T;
        const auto res = lambdaZ/(RT*ln10);
        return res;
    }

    auto Eh() const -> real
    {
        const auto T = props.temperature();
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        const auto res = ln10*RT/F*pE();
        return res;
    }

 };

DoubleLayerProps::DoubleLayerProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

DoubleLayerProps::DoubleLayerProps(const ChemicalState& state)
: pimpl(new Impl(state))
{}

DoubleLayerProps::DoubleLayerProps(ChemicalProps const& props)
: pimpl(new Impl(props))
{}

DoubleLayerProps::DoubleLayerProps(const DoubleLayerProps& other)
: pimpl(new Impl(*other.pimpl))
{}

DoubleLayerProps::~DoubleLayerProps()
{}

auto DoubleLayerProps::operator=(DoubleLayerProps other) -> DoubleLayerProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto DoubleLayerProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto DoubleLayerProps::update(ChemicalProps const& props) -> void
{
    pimpl->update(props);
}

auto DoubleLayerProps::temperature() const -> real
{
    return pimpl->temperature();
}

auto DoubleLayerProps::pressure() const -> real
{
    return pimpl->pressure();
}

auto DoubleLayerProps::elementMolality(const StringOrIndex& symbol) const -> real
{
    return pimpl->elementMolality(symbol);
}

auto DoubleLayerProps::elementMolalities() const -> ArrayXr
{
    return pimpl->elementMolalities();
}

auto DoubleLayerProps::speciesMolality(const StringOrIndex& name) const -> real
{
    return pimpl->speciesMolality(name);
}

auto DoubleLayerProps::speciesMolalities() const -> ArrayXr
{
    return pimpl->speciesMolalities();
}

auto DoubleLayerProps::speciesCharges() const -> ArrayXr
{
    return pimpl->speciesCharges();
}

auto DoubleLayerProps::ionicStrength() const -> real
{
    return pimpl->ionicStrength();
}

auto DoubleLayerProps::ionicStrengthEffective() const -> real
{
    return pimpl->ionicStrength();
}

auto DoubleLayerProps::ionicStrengthStoichiometric() const -> real
{
    return pimpl->ionicStrengthStoichiometric();
}

auto DoubleLayerProps::pH() const -> real
{
    return pimpl->pH();
}

auto DoubleLayerProps::pE() const -> real
{
    return pimpl->pE();
}

auto DoubleLayerProps::Eh() const -> real
{
    return pimpl->Eh();
}

auto DoubleLayerProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

auto DoubleLayerProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto DoubleLayerProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const DoubleLayerProps& props) -> std::ostream&
{
    const auto elements = props.phase().elements();
    const auto species = props.phase().species();
    const auto ms = props.speciesMolalities();
    const auto me = props.elementMolalities();
    const auto z = props.speciesCharges();

    assert(species.size() == ms.size());
    assert(elements.size() == me.size());
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(props.temperature()), "K" });
    table.add_row({ "Pressure", strfix(props.pressure()*1e-5), "bar" });
    table.add_row({ "Ionic Strength (Effective)", strfix(props.ionicStrength()), "molal" });
    table.add_row({ "Ionic Strength (Stoichiometric)", strfix(props.ionicStrengthStoichiometric()), "molal" });
    table.add_row({ "pH", strfix(props.pH()), "" });
    table.add_row({ "pE", strfix(props.pE()), "" });
    table.add_row({ "Eh", strfix(props.Eh()), "V" });
    table.add_row({ "Z", strfix((z*ms).sum()), "eq" });
    table.add_row({ "Element Molality:" });
    for(auto i = 0; i < elements.size(); ++i)
        if(elements[i].symbol() != "H" && elements[i].symbol() != "O")
            table.add_row({ ":: " + elements[i].symbol(), strsci(me[i]), "molal" });
    table.add_row({ "Species Molality:" });
    for(auto i = 0; i < species.size(); ++i)
        if(species[i].formula().str() != "H2O")
            table.add_row({ ":: " + species[i].repr(), strsci(ms[i]), "molal" });

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
