// Reaktoro is a unified framework for modeling chemically reactive phases.
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

#include "ChemicalProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

ChemicalProps::ChemicalProps(const ChemicalSystem& system)
: msystem(system)
{
    const auto N = system.species().size();
    const auto K = system.phases().size();

    Ts   = ArrayXr::Zero(K);
    Ps   = ArrayXr::Zero(K);
    n    = ArrayXr::Zero(N);
    nsum = ArrayXr::Zero(K);
    x    = ArrayXr::Zero(N);
    G0   = ArrayXr::Zero(N);
    H0   = ArrayXr::Zero(N);
    V0   = ArrayXr::Zero(N);
    Cp0  = ArrayXr::Zero(N);
    Cv0  = ArrayXr::Zero(N);
    Vex  = ArrayXr::Zero(K);
    VexT = ArrayXr::Zero(K);
    VexP = ArrayXr::Zero(K);
    Gex  = ArrayXr::Zero(K);
    Hex  = ArrayXr::Zero(K);
    Cpex = ArrayXr::Zero(K);
    Cvex = ArrayXr::Zero(K);
    ln_g = ArrayXr::Zero(N);
    ln_a = ArrayXr::Zero(N);
    u    = ArrayXr::Zero(N);
}

ChemicalProps::ChemicalProps(const ChemicalState& state)
: ChemicalProps(state.system())
{
    update(state);
}

auto ChemicalProps::update(const ChemicalState& state) -> void
{
    const auto& T = state.temperature();
    const auto& P = state.pressure();
    const auto& n = state.speciesAmounts();
    update(T, P, n);
}

auto ChemicalProps::update(const real& T, const real& P, ArrayXrConstRef n) -> void
{
    this->T = T;
    this->P = P;
    const auto K = msystem.phases().size();
    auto offset = 0;
    for(auto i = 0; i < K; ++i)
    {
        const auto size = msystem.phase(i).species().size();
        const auto np = n.segment(offset, size);
        phaseProps(i).update(T, P, np);
        offset += size;
    }
}

auto ChemicalProps::update(ArrayXrConstRef data) -> void
{
    ArraySerialization::deserialize(data, T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::update(ArrayXdConstRef data) -> void
{
    ArraySerialization::deserialize(data, T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::updateIdeal(const ChemicalState& state) -> void
{
    const auto& T = state.temperature();
    const auto& P = state.pressure();
    const auto& n = state.speciesAmounts();
    updateIdeal(T, P, n);
}

auto ChemicalProps::updateIdeal(const real& T, const real& P, ArrayXrConstRef n) -> void
{
    this->T = T;
    this->P = P;
    const auto K = msystem.phases().size();
    auto offset = 0;
    for(auto i = 0; i < K; ++i)
    {
        const auto size = msystem.phase(i).species().size();
        const auto np = n.segment(offset, size);
        phaseProps(i).updateIdeal(T, P, np);
        offset += size;
    }
}

auto ChemicalProps::serialize(ArrayStream<real>& stream) const -> void
{
    stream.from(T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::serialize(ArrayStream<double>& stream) const -> void
{
    stream.from(T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::deserialize(const ArrayStream<real>& stream) -> void
{
    stream.to(T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::deserialize(const ArrayStream<double>& stream) -> void
{
    stream.to(T, P, Ts, Ps, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, u);
}

auto ChemicalProps::system() const -> const ChemicalSystem&
{
    return msystem;
}

auto ChemicalProps::phaseProps(Index iphase) const -> ChemicalPropsPhaseConstRef
{
    return const_cast<ChemicalProps&>(*this).phaseProps(iphase);
}

auto ChemicalProps::phaseProps(Index iphase) -> ChemicalPropsPhaseRef
{
    const auto phase = msystem.phase(iphase);
    const auto begin = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto size = phase.species().size();

    return ChemicalPropsPhaseRef(phase, {
        Ts[iphase],
        Ps[iphase],
        n.segment(begin, size),
        nsum[iphase],
        x.segment(begin, size),
        G0.segment(begin, size),
        H0.segment(begin, size),
        V0.segment(begin, size),
        Cp0.segment(begin, size),
        Cv0.segment(begin, size),
        Vex[iphase],
        VexT[iphase],
        VexP[iphase],
        Gex[iphase],
        Hex[iphase],
        Cpex[iphase],
        Cvex[iphase],
        ln_g.segment(begin, size),
        ln_a.segment(begin, size),
        u.segment(begin, size)
    });
}

auto ChemicalProps::temperature() const -> real
{
    return T;
}

auto ChemicalProps::pressure() const -> real
{
    return P;
}

auto ChemicalProps::charge() const -> real
{
    const auto Acharge = msystem.formulaMatrixCharge();
    return Acharge * n.matrix();
}

auto ChemicalProps::elementAmount(StringOrIndex element) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto A = msystem.formulaMatrixElements();
    return A.row(ielement) * n.matrix();
}

auto ChemicalProps::elementAmountInPhase(StringOrIndex element, StringOrIndex phase) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto A = msystem.formulaMatrixElements();
    const auto np = n.matrix().segment(offset, length);
    const auto Aep = A.row(ielement).segment(offset, length);
    return Aep * np;
}

auto ChemicalProps::elementAmountAmongSpecies(StringOrIndex element, ArrayXlConstRef indices) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto A = msystem.formulaMatrixElements();
    const auto Aei = A.row(ielement)(indices);
    const auto ni = n(indices).matrix();
    return Aei * ni;
}

auto ChemicalProps::elementMass(StringOrIndex element) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmount(ielement);
    return amount * molarmass;
}

auto ChemicalProps::elementMassInPhase(StringOrIndex element, StringOrIndex phase) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmountInPhase(ielement, iphase);
    return amount * molarmass;
}

auto ChemicalProps::elementMassAmongSpecies(StringOrIndex element, ArrayXlConstRef indices) const -> real
{
    const auto ielement = detail::resolveElementIndex(msystem, element);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmountAmongSpecies(ielement, indices);
    return amount * molarmass;
}

auto ChemicalProps::speciesAmount(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return n[ispecies];
}

auto ChemicalProps::speciesMass(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return n[ispecies] * msystem.species(ispecies).molarMass();
}

auto ChemicalProps::phaseAmount(StringOrIndex phase) const -> real
{
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    return phaseProps(iphase).amount();
}

auto ChemicalProps::phaseMass(StringOrIndex phase) const -> real
{
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    return phaseProps(iphase).mass();
}

auto ChemicalProps::phaseVolume(StringOrIndex phase) const -> real
{
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    return phaseProps(iphase).volume();
}

auto ChemicalProps::moleFraction(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return x[ispecies];
}

auto ChemicalProps::concentration(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return exp(ln_a[ispecies] - ln_g[ispecies]);
}

auto ChemicalProps::lgConcentration(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return (ln_a[ispecies] - ln_g[ispecies])/ln10;
}

auto ChemicalProps::lnConcentration(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return ln_a[ispecies] - ln_g[ispecies];
}

auto ChemicalProps::activityCoefficient(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return exp(ln_g[ispecies]);
}

auto ChemicalProps::lgActivityCoefficient(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return ln_g[ispecies]/ln10;
}

auto ChemicalProps::lnActivityCoefficient(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return ln_g[ispecies];
}

auto ChemicalProps::activity(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return exp(ln_a[ispecies]);
}

auto ChemicalProps::lgActivity(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return ln_a[ispecies]/ln10;
}

auto ChemicalProps::lnActivity(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return ln_a[ispecies];
}

auto ChemicalProps::chemicalPotential(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return u[ispecies];
}

auto ChemicalProps::standardVolume(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return V0[ispecies];
}

auto ChemicalProps::standardGibbsEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return G0[ispecies];
}

auto ChemicalProps::standardEnthalpy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return H0[ispecies];
}

auto ChemicalProps::standardEntropy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return (H0[ispecies] - G0[ispecies])/T; // from G0 = H0 - T*S0
}

auto ChemicalProps::standardInternalEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return H0[ispecies] - P*V0[ispecies]; // from H0 = U0 + P*V0
}

auto ChemicalProps::standardHelmholtzEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return G0[ispecies] - P*V0[ispecies]; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ChemicalProps::standardHeatCapacityConstP(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return Cp0[ispecies];
}

auto ChemicalProps::standardHeatCapacityConstV(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndex(msystem, species);
    return Cv0[ispecies];
}

auto ChemicalProps::elementAmounts() const -> ArrayXr
{
    const auto A = msystem.formulaMatrixElements();
    return (A * n.matrix()).array();
}

auto ChemicalProps::elementAmountsInPhase(StringOrIndex phase) const -> ArrayXr
{
    const auto iphase = detail::resolvePhaseIndex(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto A = msystem.formulaMatrixElements();
    const auto Ap = A.middleCols(offset, length);
    const auto np = n.matrix().segment(offset, length);
    return (Ap * np).array();
}

auto ChemicalProps::elementAmountsAmongSpecies(ArrayXlConstRef indices) const -> ArrayXr
{
    const auto A = msystem.formulaMatrixElements();
    const auto Ai = A(Eigen::all, indices);
    const auto ni = n(indices).matrix();
    return (Ai * ni).array();
}

auto ChemicalProps::speciesAmounts() const -> ArrayXrConstRef
{
    return n;
}

auto ChemicalProps::speciesMasses() const -> ArrayXr
{
    ArrayXr m(n);
    auto i = 0; for(const auto& species : system().species())
        m[i++] *= species.molarMass();
    return m;
}

auto ChemicalProps::moleFractions() const -> ArrayXrConstRef
{
    return x;
}

auto ChemicalProps::lnConcentrations() const -> ArrayXr
{
    return ln_a - ln_g;
}

auto ChemicalProps::lnActivityCoefficients() const -> ArrayXrConstRef
{
    return ln_g;
}

auto ChemicalProps::lnActivities() const -> ArrayXrConstRef
{
    return ln_a;
}

auto ChemicalProps::chemicalPotentials() const -> ArrayXrConstRef
{
    return u;
}

auto ChemicalProps::standardVolumes() const -> ArrayXrConstRef
{
    return V0;
}

auto ChemicalProps::standardGibbsEnergies() const -> ArrayXrConstRef
{
    return G0;
}

auto ChemicalProps::standardEnthalpies() const -> ArrayXrConstRef
{
    return H0;
}

auto ChemicalProps::standardEntropies() const -> ArrayXr
{
    return (H0 - G0)/T; // from G0 = H0 - T*S0
}

auto ChemicalProps::standardInternalEnergies() const -> ArrayXr
{
    return H0 - P*V0; // from H0 = U0 + P*V0
}

auto ChemicalProps::standardHelmholtzEnergies() const -> ArrayXr
{
    return G0 - P*V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ChemicalProps::standardHeatCapacitiesConstP() const -> ArrayXrConstRef
{
    return Cp0;
}

auto ChemicalProps::standardHeatCapacitiesConstV() const -> ArrayXrConstRef
{
    return Cv0;
}

auto ChemicalProps::amount() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).amount(); });
}

auto ChemicalProps::mass() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).mass(); });
}

auto ChemicalProps::volume() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).volume(); });
}

auto ChemicalProps::gibbsEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).gibbsEnergy(); });
}

auto ChemicalProps::enthalpy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).enthalpy(); });
}

auto ChemicalProps::entropy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).entropy(); });
}

auto ChemicalProps::internalEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).internalEnergy(); });
}

auto ChemicalProps::helmholtzEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).helmholtzEnergy(); });
}

auto ChemicalProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ChemicalProps::output(const String& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

ChemicalProps::operator VectorXr() const
{
    ArrayStream<real> stream;
    serialize(stream);
    return stream.data();
}

ChemicalProps::operator VectorXd() const
{
    ArrayStream<double> stream;
    serialize(stream);
    return stream.data();
}

auto operator<<(std::ostream& out, const ChemicalProps& props) -> std::ostream&
{
    const auto species = props.system().species();
    const auto elements = props.system().elements();
    const auto b   = props.elementAmounts();
    const auto n   = props.speciesAmounts();
    const auto x   = props.moleFractions();
    const auto lng = props.lnActivityCoefficients();
    const auto lna = props.lnActivities();
    const auto mu  = props.chemicalPotentials();
    const auto G0  = props.standardGibbsEnergies();
    const auto H0  = props.standardEnthalpies();
    const auto V0  = props.standardVolumes();
    const auto S0  = props.standardEntropies();
    const auto U0  = props.standardInternalEnergies();
    const auto A0  = props.standardHelmholtzEnergies();
    const auto Cp0 = props.standardHeatCapacitiesConstP();
    const auto Cv0 = props.standardHeatCapacitiesConstV();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", str(props.temperature()), "K" });
    table.add_row({ "Pressure", str(props.pressure()), "Pa" });
    table.add_row({ "Volume", str(props.volume()), "m3" });
    table.add_row({ "Gibbs Energy", str(props.gibbsEnergy()), "J" });
    table.add_row({ "Enthalpy", str(props.enthalpy()), "J" });
    table.add_row({ "Entropy", str(props.entropy()), "J/K" });
    table.add_row({ "Internal Energy", str(props.internalEnergy()), "J" });
    table.add_row({ "Helmholtz Energy", str(props.helmholtzEnergy()), "J" });

    table.add_row({ "Element Amount:" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), str(b[i]), "mol" });
    table.add_row({ "Species Amount:" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(n[i]), "mol" });
    table.add_row({ "Mole Fraction:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(x[i]), "mol/mol" });
    table.add_row({ "Activity Coefficient:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(exp(lng[i])), "-" });
    table.add_row({ "Activity:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(exp(lna[i])), "-" });
    table.add_row({ "lg(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(lna[i]/ln10), "-" });
    table.add_row({ "ln(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(lna[i]), "-" });
    table.add_row({ "Chemical Potential:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(mu[i]), "J/mol" });
    table.add_row({ "Standard Volume:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(V0[i]), "m3/mol" });
    table.add_row({ "Standard Gibbs Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(G0[i]), "J/mol" });
    table.add_row({ "Standard Enthalpy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(H0[i]), "J/mol" });
    table.add_row({ "Standard Entropy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(S0[i]), "J/(mol*K)" });
    table.add_row({ "Standard Internal Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(U0[i]), "J/mol" });
    table.add_row({ "Standard Helmholtz Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(A0[i]), "J/mol" });
    table.add_row({ "Standard Heat Capacity (constant P):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(Cp0[i]), "J/(mol*K)" });
    table.add_row({ "Standard Heat Capacity (constant V):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(Cv0[i]), "J/(mol*K)" });

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
