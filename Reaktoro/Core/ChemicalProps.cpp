// Reaktoro is a unified framework for modeling chemically reactive phases.
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

#include "ChemicalProps.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/PhaseChemicalProps.hpp>

namespace Reaktoro {

ChemicalProps::ChemicalProps(const ChemicalSystem& system)
: sys(system)
{
    const auto numspecies = sys.species().size();
    const auto numphases = sys.phases().size();

    props.n    = ArrayXr::Zero(numspecies);
    props.nsum = ArrayXr::Zero(numphases);
    props.x    = ArrayXr::Zero(numspecies);
    props.G0   = ArrayXr::Zero(numspecies);
    props.H0   = ArrayXr::Zero(numspecies);
    props.V0   = ArrayXr::Zero(numspecies);
    props.Cp0  = ArrayXr::Zero(numspecies);
    props.Cv0  = ArrayXr::Zero(numspecies);
    props.Vex  = ArrayXr::Zero(numphases);
    props.VexT = ArrayXr::Zero(numphases);
    props.VexP = ArrayXr::Zero(numphases);
    props.Gex  = ArrayXr::Zero(numphases);
    props.Hex  = ArrayXr::Zero(numphases);
    props.Cpex = ArrayXr::Zero(numphases);
    props.Cvex = ArrayXr::Zero(numphases);
    props.ln_g = ArrayXr::Zero(numspecies);
    props.ln_a = ArrayXr::Zero(numspecies);
}

auto ChemicalProps::update(const real& T, const real& P, ArrayXrConstRef n) -> void
{
    const auto numphases = sys.phases().size();
    auto offset = 0;
    for(auto i = 0; i < numphases; ++i)
    {
        const auto size = sys.phase(i).species().size();
        const auto np = n.segment(offset, size);
        phaseProps(i).update(T, P, np);
        offset += size;
    }
}

auto ChemicalProps::update(const real& T, const real& P, ArrayXrConstRef n, Wrt<real&> wrtvar) -> void
{
    autodiff::seed(wrtvar);
    update(T, P, n);
    autodiff::unseed(wrtvar);
}

auto ChemicalProps::system() const -> const ChemicalSystem&
{
    return sys;
}

auto ChemicalProps::data() const -> const ChemicalPropsData&
{
    return props;
}

auto ChemicalProps::phaseProps(Index idx) const -> PhaseChemicalPropsConstRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.phases().numSpeciesUntilPhase(idx);
    const auto size = phase.species().size();

    return PhaseChemicalPropsConstRef(phase, {
        props.T,
        props.P,
        props.n.segment(begin, size),
        props.nsum[idx],
        props.x.segment(begin, size),
        props.G0.segment(begin, size),
        props.H0.segment(begin, size),
        props.V0.segment(begin, size),
        props.Cp0.segment(begin, size),
        props.Cv0.segment(begin, size),
        props.Vex[idx],
        props.VexT[idx],
        props.VexP[idx],
        props.Gex[idx],
        props.Hex[idx],
        props.Cpex[idx],
        props.Cvex[idx],
        props.ln_g.segment(begin, size),
        props.ln_a.segment(begin, size)
    });
}

auto ChemicalProps::phaseProps(Index idx) -> PhaseChemicalPropsRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.phases().numSpeciesUntilPhase(idx);
    const auto size = phase.species().size();

    return PhaseChemicalPropsRef(phase, {
        props.T,
        props.P,
        props.n.segment(begin, size),
        props.nsum[idx],
        props.x.segment(begin, size),
        props.G0.segment(begin, size),
        props.H0.segment(begin, size),
        props.V0.segment(begin, size),
        props.Cp0.segment(begin, size),
        props.Cv0.segment(begin, size),
        props.Vex[idx],
        props.VexT[idx],
        props.VexP[idx],
        props.Gex[idx],
        props.Hex[idx],
        props.Cpex[idx],
        props.Cvex[idx],
        props.ln_g.segment(begin, size),
        props.ln_a.segment(begin, size)
    });
}

auto ChemicalProps::temperature() const -> real
{
    return props.T;
}

auto ChemicalProps::pressure() const -> real
{
    return props.P;
}

auto ChemicalProps::speciesAmounts() const -> ArrayXrConstRef
{
    return props.n;
}

auto ChemicalProps::moleFractions() const -> ArrayXrConstRef
{
    return props.x;
}

auto ChemicalProps::lnActivityCoefficients() const -> ArrayXrConstRef
{
    return props.ln_g;
}

auto ChemicalProps::lnActivities() const -> ArrayXrConstRef
{
    return props.ln_a;
}

auto ChemicalProps::chemicalPotentials() const -> ArrayXr
{
    const auto R = universalGasConstant;
    return props.G0 + R*props.T * props.ln_a;
}

auto ChemicalProps::standardGibbsEnergies() const -> ArrayXrConstRef
{
    return props.G0;
}

auto ChemicalProps::standardEnthalpies() const -> ArrayXrConstRef
{
    return props.H0;
}

auto ChemicalProps::standardVolumes() const -> ArrayXrConstRef
{
    return props.V0;
}

auto ChemicalProps::standardEntropies() const -> ArrayXr
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

auto ChemicalProps::standardInternalEnergies() const -> ArrayXr
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

auto ChemicalProps::standardHelmholtzEnergies() const -> ArrayXr
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ChemicalProps::standardHeatCapacitiesConstP() const -> ArrayXrConstRef
{
    return props.Cp0;
}

auto ChemicalProps::standardHeatCapacitiesConstV() const -> ArrayXrConstRef
{
    return props.Cv0;
}

auto ChemicalProps::gibbsEnergy() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).gibbsEnergy(); });
}

auto ChemicalProps::enthalpy() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).enthalpy(); });
}

auto ChemicalProps::volume() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).volume(); });
}

auto ChemicalProps::entropy() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).entropy(); });
}

auto ChemicalProps::internalEnergy() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).internalEnergy(); });
}

auto ChemicalProps::helmholtzEnergy() const -> real
{
    const auto iend = sys.phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).helmholtzEnergy(); });
}

} // namespace Reaktoro
