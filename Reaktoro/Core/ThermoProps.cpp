// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright © 2014-2024 Allan Leal
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

#include "ThermoProps.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>

namespace Reaktoro {

struct ThermoProps::Impl
{
    /// The chemical system associated with these standard thermodynamic properties.
    ChemicalSystem system;

    /// The temperature of the system (in K).
    real T = 0.0;

    /// The pressure of the system (in Pa).
    real P = 0.0;

    /// The temperatures of each phase (in K).
    ArrayXr Ts;

    /// The pressures of each phase (in Pa).
    ArrayXr Ps;

    /// The standard molar Gibbs energies of formation of the species in the system (in J/mol).
    ArrayXr G0;

    /// The standard molar enthalpies of formation of the species in the system (in J/mol).
    ArrayXr H0;

    /// The standard molar volumes of the species in the system (in m3/mol).
    ArrayXr V0;

    /// The temperature derivative of the standard molar volumes of the species in the system (in m³/(mol·K)).
    ArrayXr VT0;

    /// The pressure derivative of the standard molar volumes of the species in the system (in m³/(mol·Pa)).
    ArrayXr VP0;

    /// The standard molar isobaric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cp0;

    /// Construct a ThermoProps::Impl object.
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        const auto numspecies = system.species().size();
        const auto numphases = system.phases().size();

        Ts   = ArrayXr::Zero(numphases);
        Ps   = ArrayXr::Zero(numphases);
        G0   = ArrayXr::Zero(numspecies);
        H0   = ArrayXr::Zero(numspecies);
        V0   = ArrayXr::Zero(numspecies);
        VT0  = ArrayXr::Zero(numspecies);
        VP0  = ArrayXr::Zero(numspecies);
        Cp0  = ArrayXr::Zero(numspecies);
    }

    /// Update the standard thermodynamic properties of the chemical system.
    auto update(const real& T, const real& P) -> void
    {
        this->T = T;
        this->P = P;
        const auto numphases = system.phases().size();
        for(auto i = 0; i < numphases; ++i)
            phaseProps(i).update(T, P);
    }

    /// Update the standard thermodynamic properties of the chemical system.
    auto update(const real& T, const real& P, Wrt<real&> wrtvar) -> void
    {
        autodiff::seed(wrtvar);
        update(T, P);
        autodiff::unseed(wrtvar);
    }

    /// Return the standard thermodynamic properties of a phase with given index.
    auto phaseProps(Index idx) -> ThermoPropsPhaseRef
    {
        const auto phase = system.phase(idx);
        const auto begin = system.phases().numSpeciesUntilPhase(idx);
        const auto size = phase.species().size();

        return ThermoPropsPhaseRef(phase, {
            Ts[idx],
            Ps[idx],
            G0.segment(begin, size),
            H0.segment(begin, size),
            V0.segment(begin, size),
            VT0.segment(begin, size),
            VP0.segment(begin, size),
            Cp0.segment(begin, size),
        });
    }
};

ThermoProps::ThermoProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ThermoProps::ThermoProps(const ThermoProps& other)
: pimpl(new Impl(*other.pimpl))
{}

ThermoProps::~ThermoProps()
{}

auto ThermoProps::operator=(ThermoProps other) -> ThermoProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ThermoProps::update(const real& T, const real& P) -> void
{
    pimpl->update(T, P);
}

auto ThermoProps::update(const real& T, const real& P, Wrt<real&> wrtvar) -> void
{
    pimpl->update(T, P, wrtvar);
}

auto ThermoProps::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ThermoProps::phaseProps(Index idx) const -> ThermoPropsPhaseConstRef
{
    return pimpl->phaseProps(idx);
}

auto ThermoProps::temperature() const -> const real&
{
    return pimpl->T;
}

auto ThermoProps::pressure() const -> const real&
{
    return pimpl->P;
}

auto ThermoProps::speciesStandardVolumes() const -> ArrayXrConstRef
{
    return pimpl->V0;
}

auto ThermoProps::speciesStandardVolumesT() const -> ArrayXrConstRef
{
    return pimpl->VT0;
}

auto ThermoProps::speciesStandardVolumesP() const -> ArrayXrConstRef
{
    return pimpl->VP0;
}

auto ThermoProps::speciesStandardGibbsEnergies() const -> ArrayXrConstRef
{
    return pimpl->G0;
}

auto ThermoProps::speciesStandardEnthalpies() const -> ArrayXrConstRef
{
    return pimpl->H0;
}

auto ThermoProps::speciesStandardEntropies() const -> ArrayXr
{
    return (pimpl->H0 - pimpl->G0)/pimpl->T; // from G0 = H0 - T*S0
}

auto ThermoProps::speciesStandardInternalEnergies() const -> ArrayXr
{
    return pimpl->H0 - pimpl->P * pimpl->V0; // from H0 = U0 + P*V0
}

auto ThermoProps::speciesStandardHelmholtzEnergies() const -> ArrayXr
{
    return pimpl->G0 - pimpl->P * pimpl->V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ThermoProps::speciesStandardHeatCapacitiesConstP() const -> ArrayXrConstRef
{
    return pimpl->Cp0;
}

auto ThermoProps::speciesStandardHeatCapacitiesConstV() const -> ArrayXrConstRef
{
    return pimpl->Cp0 + pimpl->T * pimpl->VT0 * pimpl->VT0 / pimpl->VP0; // from Cv0 = Cp0 + T*VT0*VT0/VP0
}

} // namespace Reaktoro
