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

#include "AqueousProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Optima includes
#include <Optima/Echelonizer.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {
namespace  {

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

    /// The formula matrix of the aqueous species.
    MatrixXd Aaq;

    /// The amounts of the species in the aqueous phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    VectorXd naq;

    /// The echelon form of the formula matrix of the aqueous species.
    Optima::Echelonizer echelonizer;

    /// The extra properties and data produced during the evaluation of the aqueous phase activity model.
    Map<String, Any> m_extra;

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

        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto A = system.formulaMatrix();

        Aaq = A.middleCols(ifirst, size);

        aqstate.T = NaN;
        aqstate.P = NaN;
        aqstate.rho = NaN;
        aqstate.epsilon = NaN;
        aqstate.Ie = NaN;
        aqstate.Is = NaN;
        aqstate.m.setConstant(size, NaN);
        aqstate.ms.setConstant(size, NaN);

        echelonizer.compute(Aaq); // echelon form of Aaq (columns of A corresponding to aqueous species)
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

    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        props.update(T, P, n.segment(ifirst, size), m_extra);
        update(props);
    }

    auto update(const ChemicalPropsPhase& aqprops) -> void
    {
        props = aqprops;
        const auto& T = props.temperature();
        const auto& P = props.pressure();
        const auto& x = props.moleFractions();
        aqstate = aqsolution.state(T, P, x);
        naq = props.speciesAmounts();
        echelonizer.updateWithPriorityWeights(naq);
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

    auto elementMolality(const String& symbol) const -> real
    {
        const auto idx = system.elements().indexWithSymbol(symbol);
        const auto& m = aqstate.m.matrix();
        return Aaq.row(idx) * m;
    }

    auto elementMolalities() const -> VectorXr
    {
        const auto E = system.elements().size();
        const auto& m = aqstate.m.matrix();
        return Aaq.topRows(E) * m;
    }

    auto speciesMolality(const String& name) const -> real
    {
        const auto idx = phase.species().indexWithName(name);
        return aqstate.m[idx];
    }

    auto speciesMolalities() const -> VectorXr
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
        const auto ln_aH = props.lnActivities()[iH];
        return -ln_aH/ln10;
    }

    auto pE() const -> real
    {
        const auto T = props.temperature();
        const auto u = props.chemicalPotentials();
        const auto ib = echelonizer.indicesBasicVariables();
        const auto R = echelonizer.R();
        const VectorXr ub = u(ib);
        const auto lambda = -R.transpose() * ub;
        const auto E = system.elements().size();
        const auto lambda_Z = lambda[E];
        const auto RT = universalGasConstant * T;
        const auto res = lambda_Z/(RT*ln10);
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

auto AqueousProps::elementMolality(const String& symbol) const -> real
{
    return pimpl->elementMolality(symbol);
}

auto AqueousProps::elementMolalities() const -> VectorXr
{
    return pimpl->elementMolalities();
}

auto AqueousProps::speciesMolality(const String& name) const -> real
{
    return pimpl->speciesMolality(name);
}

auto AqueousProps::speciesMolalities() const -> VectorXr
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
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", str(props.temperature()), "K" });
    table.add_row({ "Pressure", str(props.pressure()), "Pa" });
    table.add_row({ "Ionic Strength (Effect.)", str(props.ionicStrength()), "molal" });
    table.add_row({ "Ionic Strength (Stoich.)", str(props.ionicStrengthStoichiometric()), "molal" });
    table.add_row({ "pH", str(props.pH()), "" });
    table.add_row({ "pE", str(props.pE()), "" });
    table.add_row({ "Eh", str(props.Eh()), "V" });
    table.add_row({ "Element Molality:" });
    for(auto i = 0; i < me.size(); ++i)
        if(elements[i].symbol() != "H" && elements[i].symbol() != "O")
            table.add_row({ ":: " + elements[i].symbol(), str(me[i]), "molal" });
    table.add_row({ "Species Molality:" });
    for(auto i = 0; i < ms.size(); ++i)
        if(species[i].formula().str() != "H2O")
            table.add_row({ ":: " + species[i].name(), str(ms[i]), "molal" });

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
