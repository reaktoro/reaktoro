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
    const AqueousMixture ddlsolution;

    /// The index of the aqueous solvent species H2O in the double layer phase (not in the system!)
    const Index iH2O;

    /// The index of the aqueous solute species H+ in the double layer phase (not in the system!)
    const Index iH;

    /// The chemical properties of the system.
    ChemicalProps props;

    /// The state of the aqueous solution.
    AqueousMixtureState ddlstate;

    /// The amounts of the species in the double layer phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    VectorXd nddl;

    /// The more-fractions of the species in the double layer phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    VectorXd xddl;

    /// The chemical potentials of the elements in the double layer phase.
    VectorXr lambda;

    /// The formula matrix of the aqueous species in the double layer phase.
    MatrixXd Aaqs;

    /// The echelon form of the formula matrix `Aaqs` of the aqueous species.
    Optima::Echelonizer echelonizer;

    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexDoubleLayerPhase(system)),
      phase(system.phase(iphase)),
      props(system),
      ddlsolution(phase.species()),
      iH2O(phase.species().findWithFormula("H2O")),
      iH(phase.species().findWithFormula("H+"))
    {
        const auto size = phase.species().size();

        error(iH2O >= size, "Cannot create DoubleLayerProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H2O.");

        error(iH >= size, "Cannot create DoubleLayerProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H+.");

        // The symbols of the elements in the double layer phase
        const auto symbols = vectorize(phase.elements(), RKT_LAMBDA(x, x.symbol()));

        // Collect the species from the database that contains the elements in the double layer phase
        const auto species_same_elements = system.database().species().withElements(symbols);

        // Assemble the formula matrices of the aqueous and non-aqueous species w.r.t. elements in the double layer phase
        Aaqs = detail::assembleFormulaMatrix(phase.species(), phase.elements());

        // Initialize the aqueous state properties
        ddlstate.T = NaN;
        ddlstate.P = NaN;
        ddlstate.rho = NaN;
        ddlstate.epsilon = NaN;
        ddlstate.Ie = NaN;
        ddlstate.Is = NaN;
        ddlstate.m.setConstant(size, NaN);
        ddlstate.ms.setConstant(size, NaN);
        ddlstate.z.setConstant(size, NaN);

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
        auto const& ddlprops = cprops.phaseProps(iphase);
        auto const& T = ddlprops.temperature();
        auto const& P = ddlprops.pressure();
        auto const& x = ddlprops.speciesMoleFractions();

        // Update the internal properties of the double layer phase
        props = cprops;

        // Update the internal aqueous state object
        ddlstate = ddlsolution.state(T, P, x);

        // Update auxiliary vector nddl to be used in the echelonization below
        nddl = ddlprops.speciesAmounts();
        xddl = ddlprops.speciesAmounts();

//        std::cout << "x = " << x.transpose() << std::endl;
//        std::cout << "nddl = " << nddl.transpose() << std::endl;
//
//        std::cout << "DDL state:" << std::endl;
//        std::cout << "z = " << ddlstate.z.transpose() << std::endl;
//        std::cout << "m = " << ddlstate.m.transpose() << std::endl;
//        getchar();

        // Update the echelon form and also the list of basic species
        echelonizer.updateWithPriorityWeights(nddl);

        // Compute chemical potentials of the elements in the double layer phase
        const auto u = ddlprops.speciesChemicalPotentials();
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
        const auto& m = ddlstate.m.matrix();
        return Aaqs.row(idx) * m;
    }

    auto elementMolalities() const -> ArrayXr
    {
        const auto E = phase.elements().size();
        const auto& m = ddlstate.m.matrix();
        return Aaqs.topRows(E) * m;
    }

    auto speciesMolality(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return ddlstate.m[idx];
    }

    auto speciesMolalities() const -> ArrayXr
    {
        return ddlstate.m;
    }

    auto speciesCharges() const -> ArrayXr
    {
        return ddlstate.z;
    }

    auto speciesAmounts() const -> ArrayXr
    {
        return nddl;
    }

    auto speciesAmount(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return nddl[idx];
    }

    auto speciesMoleFractions() const -> ArrayXr
    {
        return xddl;
    }

    auto speciesMoleFraction(const StringOrIndex& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndex(phase, name);
        return xddl[idx];
    }

    auto ionicStrength() const -> real
    {
        return ddlstate.Ie;
    }

    auto ionicStrengthStoichiometric() const -> real
    {
        return ddlstate.Is;
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

auto DoubleLayerProps::speciesAmounts() const -> ArrayXr
{
    return pimpl->speciesAmounts();
}

auto DoubleLayerProps::speciesAmount(const StringOrIndex& symbol) const -> real
{
    return pimpl->speciesAmount(symbol);
}

auto DoubleLayerProps::speciesMoleFractions() const -> ArrayXr
{
    return pimpl->speciesMoleFractions();
}

auto DoubleLayerProps::speciesMoleFraction(const StringOrIndex& symbol) const -> real
{
    return pimpl->speciesAmount(symbol);
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
    const auto ns = props.speciesAmounts();
    const auto xs = props.speciesMoleFractions();

    assert(species.size() == ms.size());
    assert(elements.size() == me.size());
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(props.temperature()), "K" });
    table.add_row({ "Pressure", strfix(props.pressure()*1e-5), "bar" });
    table.add_row({ "Species Amounts:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), strsci(ns[i]), "moles" });
    table.add_row({ "Species Mole Fractions:" });
    for(auto i = 0; i < species.size(); ++i)
        table.add_row({ ":: " + species[i].repr(), strsci(xs[i]), "" });

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
