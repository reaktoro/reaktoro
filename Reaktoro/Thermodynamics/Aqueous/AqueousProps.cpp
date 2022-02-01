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

#include "AqueousProps.hpp"

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

/// Return the index of the first aqueous phase in the system.
auto indexAqueousPhase(const ChemicalSystem& system) -> Index
{
    const auto aqueous_phases = system.phases().withAggregateState(AggregateState::Aqueous);
    warning(aqueous_phases.size() > 1,
        "While creating an AqueousProps object, it has been detected ",
        "more than one aqueous phase in the system. The AqueousProps object "
        "created will correspond to the first aqueous phase found.");
    const auto idx = system.phases().findWithAggregateState(AggregateState::Aqueous);
    error(idx >= system.phases().size(),
        "Could not create an AqueousProps object because there is no "
        "phase in the system with aggregate state value AggregateState::Aqueous.");
    return idx;
}

// Return a chemical potential function for a species using a given activity model.
auto chemicalPotentialModel(Species const& species, ActivityModelGenerator const& generator) -> Fn<real(real, real)>
{
    const auto activitymodel = generator({species}); // TODO: Use .withMemoization() here to avoid full recomputations when same T and P are given.
    const auto R = universalGasConstant;
    const auto x = ArrayXr{{1.0}}; // the mole fraction of the single species in a pure phase
    auto actprops = ActivityProps::create(1);

    return [=](real T, real P) mutable -> real
    {
        activitymodel(actprops, {T, P, x}); // evaluate the activity model
        const auto G0 = species.standardThermoProps(T, P).G0;
        const auto ln_a = actprops.ln_a[0];
        return G0 + R*T*ln_a;
    };
}

// Return a default chemical potential function for a fluid species using ideal
// activity model. This chemical potential model assumes the species constitute
// a pure ideal gas phase (containing just a single gaseous species).
auto defaultChemicalPotentialModelFluidSpecies(Species const& species) -> Fn<real(real, real)>
{
    return chemicalPotentialModel(species, ActivityModelIdealGas());
}

// Return a default chemical potential function for a solid species using ideal
// activity model for solid solutions. This chemical potential model assumes
// the species constitute a pure solid phase (containing just a single solid
// species). In this case, an ideal model coincides with the fact that the
// activity of a pure solid is one and its chemical potential is identical to
// its standard chemical potential.
auto defaultChemicalPotentialModelSolidSpecies(Species const& species) -> Fn<real(real, real)>
{
    return [=](real T, real P)
    {
        return species.standardThermoProps(T, P).G0;
    };
}

// Return a default chemical potential function for a chemical species.
auto defaultChemicalPotentialModel(Species const& species) -> Fn<real(real, real)>
{
    switch(species.aggregateState())
    {
    case AggregateState::Gas:
    case AggregateState::Liquid:
    case AggregateState::Fluid:
        return defaultChemicalPotentialModelFluidSpecies(species);
    default:
        return defaultChemicalPotentialModelSolidSpecies(species);
    }
}

// Return a vector with default chemical potential functions for every given chemical species.
auto defaultChemicalPotentialModels(SpeciesList const& nonaqueous) -> Vec<Fn<real(real, real)>>
{
    return vectorize(nonaqueous, RKT_LAMBDA(x, defaultChemicalPotentialModel(x)));
}

} // namespace

struct AqueousProps::Impl
{
    /// The chemical system in which the aqueous phase is.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the aqueous phase in the system.
    const Index iphase;

    /// The underlying Phase object for the aqueous phase in the system.
    const Phase phase;

    /// The phase as an aqueous solution.
    const AqueousMixture aqsolution;

    /// The index of the aqueous solvent species H2O
    const Index iH2O;

    /// The index of the aqueous solute species H+
    const Index iH;

    /// The chemical properties of the aqueous phase.
    ChemicalPropsPhase props;

    /// The state of the aqueous solution.
    AqueousMixtureState aqstate;

    /// The amounts of the species in the aqueous phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
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
    Vec<Fn<real(real, real)>> chemical_potential_models;

    /// The extra properties and data produced during the evaluation of the aqueous phase activity model.
    Map<String, Any> extra;

    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexAqueousPhase(system)),
      phase(system.phase(iphase)),
      props(phase),
      aqsolution(phase.species()),
      iH2O(phase.species().findWithFormula("H2O")),
      iH(phase.species().findWithFormula("H+"))
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == aqsolution.species().size());

        const auto size = phase.species().size();

        error(iH2O >= size, "Cannot create AqueousProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H2O.");

        error(iH >= size, "Cannot create AqueousProps object for phase ", phase.name(), " "
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

        // Initialize the chemical potential models for the non-aqueous species, as if they were pure phases
        chemical_potential_models = defaultChemicalPotentialModels(nonaqueous);

        // Initialize the aqueous state properties
        aqstate.T = NaN;
        aqstate.P = NaN;
        aqstate.rho = NaN;
        aqstate.epsilon = NaN;
        aqstate.Ie = NaN;
        aqstate.Is = NaN;
        aqstate.m.setConstant(size, NaN);
        aqstate.ms.setConstant(size, NaN);

        // Compute the initial echelon form of formula matrix `Aaqs`
        echelonizer.compute(Aaqs);
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

    auto setActivityModel(const StringOrIndex& species, const ActivityModelGenerator& generator) -> void
    {
        const auto i = detail::resolveSpeciesIndex(nonaqueous, species);
        errorif(i >= nonaqueous.size(), "It was not possible to set the activity model "
            "of species with name or index `", detail::stringfy(species), "` because "
            "there is no such species in the list of species returned by method "
            "AqueousProps::saturationSpecies.");
        chemical_potential_models[i] = chemicalPotentialModel(nonaqueous[i], generator);
    }

    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        props.update(T, P, n.segment(ifirst, size), extra);
        update(props);
    }

    auto update(const ChemicalPropsPhase& aqprops) -> void
    {
        // Auxiliary variables
        const auto& T = aqprops.temperature();
        const auto& P = aqprops.pressure();
        const auto& x = aqprops.speciesMoleFractions();

        // Update the internal properties of the aqueous phase
        props = aqprops;

        // Update the internal aqueous state object
        aqstate = aqsolution.state(T, P, x);

        // Update auxiliary vector naq to be used in the echelonization below
        naq = props.speciesAmounts();

        // Update the echelon form and also the list of basic species
        echelonizer.updateWithPriorityWeights(naq);

        // Compute chemical potentials of the elements in the aqueous phase
        const auto u = props.speciesChemicalPotentials();
        const auto ib = echelonizer.indicesBasicVariables();
        const auto R = echelonizer.R();
        const auto Rb = R.topRows(ib.size());
        const VectorXr ub = u(ib);
        lambda = Rb.transpose() * ub;
    }

    auto update(const ChemicalProps& sysprops) -> void
    {
        update(sysprops.phaseProps(iphase));
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
        const auto ln_aH = props.speciesActivitiesLn()[iH];
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

    auto alkalinity() const -> real
    {
        error(true, "AqueousProps::alkalinity has not been implemented yet.");
        return {};
    }

    auto saturationIndexLn(const StringOrIndex& species) const -> real
    {
        const auto i = detail::resolveSpeciesIndex(nonaqueous, species);
        errorif(i >= nonaqueous.size(), "It was not possible to calculate the "
            "saturation index of species with name or index `", detail::stringfy(species), "` "
            "because there is no such species in the list of species returned by method "
            "AqueousProps::saturationSpecies.");
        const auto T = temperature();
        const auto P = pressure();
        const auto RT = universalGasConstant * T;
        const auto ui = chemical_potential_models[i](T, P);
        const auto li = Anon.col(i).dot(lambda);
        const auto lnOmegai = (li - ui)/RT;
        return lnOmegai;
    }

    auto saturationIndicesLn() const -> ArrayXr
    {
        const auto T = temperature();
        const auto P = pressure();
        const auto RT = universalGasConstant * T;
        const auto num_nonaqueous = nonaqueous.size();
        ArrayXr lnOmega(num_nonaqueous);
        lnOmega = Anon.transpose() * lambda;
        for(auto i = 0; i < num_nonaqueous; ++i)
            lnOmega[i] -= chemical_potential_models[i](T, P);
        lnOmega /= RT;
        return lnOmega;
    }
};

AqueousProps::AqueousProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

AqueousProps::AqueousProps(const ChemicalState& state)
: pimpl(new Impl(state.system(), state))
{}

AqueousProps::AqueousProps(const ChemicalProps& props)
: pimpl(new Impl(props.system(), props))
{}

AqueousProps::AqueousProps(const AqueousProps& other)
: pimpl(new Impl(*other.pimpl))
{}

AqueousProps::~AqueousProps()
{}

auto AqueousProps::operator=(AqueousProps other) -> AqueousProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousProps::setActivityModel(const StringOrIndex& species, const ActivityModelGenerator& generator) -> void
{
    pimpl->setActivityModel(species, generator);
}

auto AqueousProps::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto AqueousProps::update(const ChemicalProps& props) -> void
{
    pimpl->update(props);
}

auto AqueousProps::temperature() const -> real
{
    return pimpl->temperature();
}

auto AqueousProps::pressure() const -> real
{
    return pimpl->pressure();
}

auto AqueousProps::elementMolality(const StringOrIndex& symbol) const -> real
{
    return pimpl->elementMolality(symbol);
}

auto AqueousProps::elementMolalities() const -> ArrayXr
{
    return pimpl->elementMolalities();
}

auto AqueousProps::speciesMolality(const StringOrIndex& name) const -> real
{
    return pimpl->speciesMolality(name);
}

auto AqueousProps::speciesMolalities() const -> ArrayXr
{
    return pimpl->speciesMolalities();
}

auto AqueousProps::ionicStrength() const -> real
{
    return pimpl->ionicStrength();
}

auto AqueousProps::ionicStrengthEffective() const -> real
{
    return pimpl->ionicStrength();
}

auto AqueousProps::ionicStrengthStoichiometric() const -> real
{
    return pimpl->ionicStrengthStoichiometric();
}

auto AqueousProps::pH() const -> real
{
    return pimpl->pH();
}

auto AqueousProps::pE() const -> real
{
    return pimpl->pE();
}

auto AqueousProps::Eh() const -> real
{
    return pimpl->Eh();
}

auto AqueousProps::alkalinity() const -> real
{
    return pimpl->alkalinity();
}

auto AqueousProps::saturationSpecies() const -> SpeciesList
{
    return pimpl->nonaqueous;
}

auto AqueousProps::saturationIndex(const StringOrIndex& species) const -> real
{
    return exp(pimpl->saturationIndexLn(species));
}

auto AqueousProps::saturationIndexLn(const StringOrIndex& species) const -> real
{
    return pimpl->saturationIndexLn(species);
}

auto AqueousProps::saturationIndexLg(const StringOrIndex& species) const -> real
{
    return pimpl->saturationIndexLn(species) / ln10;
}

auto AqueousProps::saturationIndices() const -> ArrayXr
{
    return saturationIndicesLn().exp();
}

auto AqueousProps::saturationIndicesLn() const -> ArrayXr
{
    return pimpl->saturationIndicesLn();
}

auto AqueousProps::saturationIndicesLg() const -> ArrayXr
{
    return saturationIndicesLn() / ln10;
}

auto AqueousProps::phase() const -> const Phase&
{
    return pimpl->phase;
}

auto AqueousProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto AqueousProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const AqueousProps& props) -> std::ostream&
{
    const auto elements = props.phase().elements();
    const auto species = props.phase().species();
    const auto ms = props.speciesMolalities();
    const auto me = props.elementMolalities();
    const auto lgOmega = props.saturationIndicesLg();
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
    table.add_row({ "Element Molality:" });
    for(auto i = 0; i < elements.size(); ++i)
        if(elements[i].symbol() != "H" && elements[i].symbol() != "O")
            table.add_row({ ":: " + elements[i].symbol(), strsci(me[i]), "molal" });
    table.add_row({ "Species Molality:" });
    for(auto i = 0; i < species.size(); ++i)
        if(species[i].formula().str() != "H2O")
            table.add_row({ ":: " + species[i].repr(), strsci(ms[i]), "molal" });
    table.add_row({ "Saturation Indices (log base 10):" });
    for(auto [i, species] : enumerate(props.saturationSpecies()))
        table.add_row({ ":: " + species.repr(), strfix(((lgOmega[i] + 1000)) - 1000), "-" }); // + 1000 - 1000 as a trick to transform -1e15 into 0.0

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
