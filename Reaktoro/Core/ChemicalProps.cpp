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
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>

namespace Reaktoro {

ChemicalProps::ChemicalProps(const ChemicalSystem& system)
: sys(system)
{
    const auto numspecies = sys.species().size();
    const auto numphases = sys.phases().size();

    props.n      = ArrayXr::Zero(numspecies);
    props.amount = ArrayXr::Zero(numphases);
    props.x      = ArrayXr::Zero(numspecies);
    props.G0     = ArrayXr::Zero(numspecies);
    props.H0     = ArrayXr::Zero(numspecies);
    props.V0     = ArrayXr::Zero(numspecies);
    props.Cp0    = ArrayXr::Zero(numspecies);
    props.Cv0    = ArrayXr::Zero(numspecies);
    props.Vex    = ArrayXr::Zero(numphases);
    props.VexT   = ArrayXr::Zero(numphases);
    props.VexP   = ArrayXr::Zero(numphases);
    props.Gex    = ArrayXr::Zero(numphases);
    props.Hex    = ArrayXr::Zero(numphases);
    props.Cpex   = ArrayXr::Zero(numphases);
    props.Cvex   = ArrayXr::Zero(numphases);
    props.ln_g   = ArrayXr::Zero(numspecies);
    props.ln_a   = ArrayXr::Zero(numspecies);
}

ChemicalProps::ChemicalProps(const ChemicalSystem& system, const ChemicalPropsData& data)
: sys(system), props(data)
{}

auto ChemicalProps::update(double T, double P, VectorXrConstRef n) -> void
{
    *this = sys.props(T, P, n);
}

auto ChemicalProps::system() const -> const ChemicalSystem&
{    return sys;
}

auto ChemicalProps::data() const -> const ChemicalPropsData&
{    return props;
}

auto ChemicalProps::phase(Index idx) const -> ChemicalPropsPhaseConstRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.indexFirstSpeciesInPhase(idx);
    const auto size = phase.species().size();

    return ChemicalPropsPhaseConstRef(phase, {
        props.T,
        props.P,
        props.n.segment(begin, size),
        props.amount[idx],
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

auto ChemicalProps::phase(Index idx) -> ChemicalPropsPhaseRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.indexFirstSpeciesInPhase(idx);
    const auto size = phase.species().size();

    return ChemicalPropsPhaseRef(phase, {
        props.T,
        props.P,
        props.n.segment(begin, size),
        props.amount[idx],
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

auto ChemicalProps::fluidVolume() const -> real
{
    real sum = 0.0;
    for(auto i = 0; i < sys.phases().size(); ++i)
        if(sys.phase(i).stateOfMatter() != StateOfMatter::Solid)
            sum += phase(i).volume();
    return sum;
}

auto ChemicalProps::solidVolume() const -> real
{
    real sum = 0.0;
    for(auto i = 0; i < sys.phases().size(); ++i)
        if(sys.phase(i).stateOfMatter() == StateOfMatter::Solid)
            sum += phase(i).volume();
    return sum;
}

auto ChemicalProps::volume() const -> real
{
    real sum = 0.0;
    for(auto i = 0; i < sys.phases().size(); ++i)
        sum += phase(i).volume();
    return sum;
}

} // namespace Reaktoro
