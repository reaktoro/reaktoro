// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>

namespace Reaktoro {
namespace detail {

/// Return the molar masses of the species
auto molarMasses(const SpeciesList& species)
{
    ArrayXd molar_masses(species.size());
    transform(species, molar_masses, [](auto&& s) { return s.molarMass(); });
    return molar_masses;
}

} // namespace detail

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The state of matter of the phase.
    StateOfMatter state = StateOfMatter::Solid;

    /// The list of Species instances defining the phase
    SpeciesList species;

    /// The standard thermodynamic model of the species in the phase.
    StandardThermoModelFn standard_thermo_model_fn;

    /// The activity model of the phase.
    ActivityPropsFn activity_props_fn;

    /// The molar masses of the species in the phase.
    ArrayXd species_molar_masses;

    /// Evaluate the standard thermodynamic properties of the phase.
    auto eval(ThermoPropsPhaseRef props, real T, real P) const -> void
    {
        auto&& data = props.data();

        const auto updatingT = data.T != T;
        const auto updatingP = data.P != P;

        if(updatingT || updatingP)
        {
            data.T = T;
            data.P = P;
            evalStandardThermoProps(props);
        }
    }

    /// Evaluate the chemical properties of the phase.
    auto eval(ChemicalPropsPhaseRef props, real T, real P, ArrayXrConstRef n) const -> void
    {
        auto&& data = props.data();

        const auto updatingT = data.T != T;
        const auto updatingP = data.P != P;
        const auto updatingN = (data.n != n).all();

        if(updatingT || updatingP)
        {
            data.T = T;
            data.P = P;
            evalStandardThermoProps(props);
        }
        if(updatingT || updatingP || updatingN)
        {
            data.n = n;
            evalMoleFractions(props, n);
            evalActivityProps(props);
        }
    }

    /// Evaluate the mole fractions of the species in the phase.
    auto evalMoleFractions(ChemicalPropsPhaseRef props, ArrayXrConstRef n) const -> void
    {
        auto&& data = props.data();
        auto& amount = data.amount;
        auto& x = data.x;
        amount = n.sum();
        const auto size = species.size();
        assert(x.size() == size);
        if(amount == 0.0)
            x = (size == 1) ? 1.0 : 0.0;
        else x = n/amount;
    }

    /// Evaluate the activity properties of the phase.
    auto evalActivityProps(ChemicalPropsPhaseRef props) const -> void
    {
        auto&& data  = props.data();

        auto& T      = data.T;
        auto& P      = data.P;
        auto& amount = data.amount;
        auto& x      = data.x;
        auto& Vex    = data.Vex;
        auto& VexT   = data.VexT;
        auto& VexP   = data.VexP;
        auto& Gex    = data.Gex;
        auto& Hex    = data.Hex;
        auto& Cpex   = data.Cpex;
        auto& Cvex   = data.Cvex;
        auto& ln_g   = data.ln_g;
        auto& ln_a   = data.ln_a;

        const auto size = species.size();
        assert(ln_g.size() == size);
        assert(ln_a.size() == size);

        if(amount == 0.0)
        {
            Vex = VexT = VexP = Gex = Hex = Cpex = Cvex = 0.0;
            ln_g = ln_a = 0.0;
        }
        else
        {
            activity_props_fn({Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a }, T, P, x);
        }
    }

    /// Evaluate the standard thermodynamic properties of the species in the phase.
    template<typename ThermoPropsRefType>
    auto evalStandardThermoProps(ThermoPropsRefType props) const -> void
    {
        auto&& data  = props.data();
        auto& T      = data.T;
        auto& P      = data.P;
        auto& G0     = data.G0;
        auto& H0     = data.H0;
        auto& V0     = data.V0;
        auto& Cp0    = data.Cp0;
        auto& Cv0    = data.Cv0;

        const auto size = species.size();
        assert(G0.size() == size);
        assert(H0.size() == size);
        assert(V0.size() == size);
        assert(Cp0.size() == size);
        assert(Cv0.size() == size);

        StandardThermoProps stdprops;
        for(auto i = 0; i < size; ++i)
        {
            stdprops = standard_thermo_model_fn(T, P, species[i]);
            G0[i]  = stdprops.G0;
            H0[i]  = stdprops.H0;
            V0[i]  = stdprops.V0;
            Cp0[i] = stdprops.Cp0;
            Cv0[i] = stdprops.Cv0;
        }
    }
};

Phase::Phase()
: pimpl(new Impl())
{}

auto Phase::withName(std::string name) -> Phase
{
    Phase copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Phase::withSpecies(SpeciesList species) -> Phase
{
    Phase copy = clone();
    copy.pimpl->species = std::move(species);
    copy.pimpl->species_molar_masses = detail::molarMasses(copy.pimpl->species);
    return copy;
}

auto Phase::withStateOfMatter(StateOfMatter state) -> Phase
{
    Phase copy = clone();
    copy.pimpl->state = std::move(state);
    return copy;
}

auto Phase::withStandardThermoModel(StandardThermoModelFn fn) -> Phase
{
    Phase copy = clone();
    copy.pimpl->standard_thermo_model_fn = std::move(fn);
    return copy;
}

auto Phase::withActivityModel(ActivityPropsFn fn) -> Phase
{
    Phase copy = clone();
    copy.pimpl->activity_props_fn = std::move(fn);
    return copy;
}

auto Phase::name() const -> std::string
{
    return pimpl->name;
}

auto Phase::stateOfMatter() const -> StateOfMatter
{
    return pimpl->state;
}

auto Phase::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto Phase::species(Index idx) const -> const Species&
{
    return pimpl->species[idx];
}

auto Phase::standardThermoModel() const -> const StandardThermoModelFn&
{
    return pimpl->standard_thermo_model_fn;
}

auto Phase::activityModel() const -> const ActivityPropsFn&
{
    return pimpl->activity_props_fn;
}

auto Phase::props(real T, real P) const -> ThermoPropsPhase
{
    ThermoPropsPhase res(*this);
    eval(res, T, P);
    return res;
}

auto Phase::props(real T, real P, ArrayXrConstRef n) const -> ChemicalPropsPhase
{
    ChemicalPropsPhase res(*this);
    eval(res, T, P, n);
    return res;
}

auto Phase::eval(ThermoPropsPhaseRef props, real T, real P) const -> void
{
    pimpl->eval(props, T, P);
}

auto Phase::eval(ChemicalPropsPhaseRef props, real T, real P, ArrayXrConstRef n) const -> void
{
    pimpl->eval(props, T, P, n);
}

auto Phase::clone() const -> Phase
{
    Phase phase;
    *phase.pimpl = *pimpl;
    return phase;
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
