// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

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
    ChemicalSystem system;

    /// The index of the underlying Phase object for the aqueous phase in the system.
    Index iphase;

    /// The underlying Phase object for the aqueous phase in the system.
    Phase phase;

    /// The chemical properties of the aqueous phase.
    ChemicalPropsPhase props;

    /// The phase as an aqueous solution.
    AqueousMixture aqsolution;

    /// The state of the aqueous solution.
    AqueousMixtureState aqstate;

    /// The index of the aqueous solvent species H2O
    Index iH2O;

    /// The index of the aqueous solute species H+
    Index iH;

    /// The formula matrix of the aqueous species.
    MatrixXd Aaq;

    /// The echelon form of the formula matrix of the aqueous species.
    Optima::Echelonizer echelonizer;

    /// Construct an AqueousProps::Impl object.
    Impl(const ChemicalSystem& system)
    : system(system),
      iphase(indexAqueousPhase(system)),
      phase(system.phase(iphase)),
      props(phase),
      aqsolution(phase.species())
    {
        assert(phase.species().size() > 0);
        assert(phase.species().size() == props.phase().species().size());
        assert(phase.species().size() == aqsolution.species().size());

        iH2O = phase.species().findWithFormula("H2O");
        iH = phase.species().findWithFormula("H+");

        const auto size = phase.species().size();

        error(iH2O >= size, "Cannot create AqueousProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H2O.");

        error(iH >= size, "Cannot create AqueousProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H+.");

        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto A = system.formulaMatrix();

        Aaq = A.middleCols(ifirst, size);

        echelonizer.compute(Aaq); // echelon form of Aaq (columns of A corresponding to aqueous species)
    }

    /// Construct an AqueousProps::Impl object.
    Impl(const ChemicalSystem& system, const ChemicalState& state)
    : Impl(system)
    {
        update(state);
    }

    /// Construct an AqueousProps::Impl object.
    Impl(const ChemicalSystem& system, const ChemicalProps& props)
    : Impl(system)
    {
        update(props);
    }

    /// Update the aqueous properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void
    {
        const auto T = state.temperature();
        const auto P = state.pressure();
        const auto n = state.speciesAmounts();
        const auto ifirst = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = phase.species().size();
        const auto naq = n.segment(ifirst, size);
        props.update(T, P, naq);
    }

    /// Update the aqueous properties with given chemical properties of the system.
    auto update(const ChemicalProps& sysprops) -> void
    {
        props.update(sysprops.phaseProps(iphase).data());
    }

    /// Return the molality of an element (in molal).
    auto elementMolality(const String& symbol) const -> real
    {
        const auto idx = system.elements().indexWithSymbol(symbol);
        const auto naq = props.speciesAmounts();
        real bi = {};
        for(auto i = 0; i < naq.size(); ++i)
            bi += Aaq(idx, i) * naq[i];
        const auto nH2O = naq[iH2O];
        const auto molality = bi/(waterMolarMass * nH2O);
        return molality;
    }

    /// Return the molality of an aqueous solute species (in molal).
    auto speciesMolality(const String& name) const -> real
    {
        const auto idx = phase.species().indexWithName(name);
        const auto naq = props.speciesAmounts();
        const auto ni = naq[idx];
        const auto nH2O = naq[iH2O];
        const auto molality = ni/(waterMolarMass * nH2O);
        return molality;
    }

    /// Return the effective ionic strength of the aqueous phase (in molal).
    auto ionicStrength() const -> real
    {
        return aqstate.Ie;
    }

    /// Return the stoichiometric ionic strength of the aqueous phase (in molal).
    auto ionicStrengthStoichiometric() const -> real
    {
        return aqstate.Is;
    }

    /// Return the pH of the aqueous phase.
    auto pH() const -> real
    {
        return props.lnActivities()[iH];
    }

    /// Return the pE of the aqueous phase.
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

    /// Return the reduction potential of the aqueous phase (in V).
    auto Eh() const -> real
    {
        const auto T = props.temperature();
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        const auto res = ln10*RT/F*pE();
        return res;
    }

    /// Return the total alkalinity of the aqueous phase (in eq/L).
    auto alkalinity() const -> real
    {
        error(true, "AqueousProps::alkalinity is has not been implemented yet.");
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

auto AqueousProps::elementMolality(const String& symbol) const -> real
{
    return pimpl->elementMolality(symbol);
}

auto AqueousProps::speciesMolality(const String& name) const -> real
{
    return pimpl->speciesMolality(name);
}

auto AqueousProps::ionicStrength() const -> real
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

} // namespace Reaktoro
