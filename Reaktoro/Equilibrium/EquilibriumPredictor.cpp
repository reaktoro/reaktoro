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

#include "EquilibriumPredictor.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>

namespace Reaktoro {
namespace {

/// The auxiliary return type for the temperature and pressure getter functions below.
using GetterFn = Fn<double(VectorXdConstRef, VectorXdConstRef)>;

/// Create the lambda function that extracts temperature from either `p` or `w`.
auto getTemperatureFn(Strings const& inputs) -> GetterFn
{
    const auto idxT = index(inputs, "T");
    const auto knownT = idxT < inputs.size(); // T is known if it is an input variable (it then lives in w)
    if(knownT) return [idxT](VectorXdConstRef p, VectorXdConstRef w) { return w[idxT]; };
    else return [](VectorXdConstRef p, VectorXdConstRef w) { return p[0]; };
}

/// Create the lambda function that extracts pressure from either `p` or `w`.
auto getPressureFn(Strings const& inputs) -> GetterFn
{
    const auto idxP = index(inputs, "P");
    const auto knownP = idxP < inputs.size(); // P is known if it is an input variable (it then lives in w)
    if(knownP) return [idxP](VectorXdConstRef p, VectorXdConstRef w) { return w[idxP]; };
    const auto idxT = index(inputs, "T");
    const auto knownT = idxT < inputs.size(); // T is known if it is an input variable (it then lives in w)
    if(knownT) return [](VectorXdConstRef p, VectorXdConstRef w) { return p[0]; }; // if T is known, then P is in p[0]
    else return [](VectorXdConstRef p, VectorXdConstRef w) { return p[1]; }; // if T is unknown, then P is in p[1]
}

} // namespace

struct EquilibriumPredictor::Impl
{
    const ChemicalState state0;
    const EquilibriumSensitivity sensitivity0;
    const VectorXd n0;    ///< The species amounts *n* at the reference equilibrium state.
    const VectorXd p0;    ///< The control variables *p* at the reference equilibrium state.
    const VectorXd q0;    ///< The control variables *q* at the reference equilibrium state.
    const VectorXd w0;    ///< The input variables *w* at the reference equilibrium state.
    const VectorXd c0;    ///< The component amounts *c* at the reference equilibrium state.
    const VectorXd u0;    ///< The chemical properties *u* at the reference equilibrium state.
    const Index Nn;       ///< The size of vector *n* with amounts of the species in the chemical system.
    const Index Nu;       ///< The size of vector *u* with the serialized properties of the chemical system.
    GetterFn getT;        ///< The function that gets temperature from either *p* or *w* depending if it is known or unwknon in the equilibrium calculation.
    GetterFn getP;        ///< The function that gets pressure from either *p* or *w* depending if it is known or unwknon in the equilibrium calculation.

    /// Construct a EquilibriumPredictor object.
    Impl(ChemicalState const& state0, EquilibriumSensitivity const& sensitivity0)
    : state0(state0), sensitivity0(sensitivity0),
      n0(state0.speciesAmounts()),
      p0(state0.equilibrium().p()),
      q0(state0.equilibrium().q()),
      w0(state0.equilibrium().w()),
      c0(state0.equilibrium().c()),
      u0(state0.props()),
      Nn(n0.size()),
      Nu(u0.size()),
      getT(getTemperatureFn(state0.equilibrium().inputNames())),
      getP(getPressureFn(state0.equilibrium().inputNames()))
    {
        errorif(state0.equilibrium().inputNames().empty(),
            "EquilibriumPredictor expects a ChemicalState object that "
            "has been used in a call to EquilibriumSolver::solve.");
    }

    auto predict(ChemicalState& state, EquilibriumConditions const& conditions) const -> void
    {
        const auto wvals = conditions.inputValues();
        const auto cvals = conditions.initialComponentAmountsGetOrCompute(state);

        const auto w = wvals.cast<double>().matrix();
        const auto c = cvals.cast<double>().matrix();

        const VectorXd dw = w - w0;
        const VectorXd dc = c - c0;

        predict(state, dw, dc);
    }

    auto predict(ChemicalState& state, VectorXdConstRef const& dw, VectorXdConstRef const& dc) const -> void
    {
        const auto dndw0 = sensitivity0.dndw(); // The derivatives *dn/dw* at the reference equilibrium state.
        const auto dpdw0 = sensitivity0.dpdw(); // The derivatives *dp/dw* at the reference equilibrium state.
        const auto dqdw0 = sensitivity0.dqdw(); // The derivatives *dq/dw* at the reference equilibrium state.
        const auto dudw0 = sensitivity0.dudw(); // The derivatives *du/dw* at the reference equilibrium state.
        const auto dndc0 = sensitivity0.dndc(); // The derivatives *dn/dc* at the reference equilibrium state.
        const auto dpdc0 = sensitivity0.dpdc(); // The derivatives *dp/dc* at the reference equilibrium state.
        const auto dqdc0 = sensitivity0.dqdc(); // The derivatives *dq/dc* at the reference equilibrium state.
        const auto dudc0 = sensitivity0.dudc(); // The derivatives *du/dc* at the reference equilibrium state.

        const auto n = n0 + dndw0*dw + dndc0*dc;
        const auto p = p0 + dpdw0*dw + dpdc0*dc;
        const auto q = q0 + dqdw0*dw + dqdc0*dc;
        const auto u = u0 + dudw0*dw + dudc0*dc;

        const auto w = w0 + dw;
        const auto c = c0 + dc;

        state.setSpeciesAmounts(n);
        state.props().update(u);
        state.equilibrium() = state0.equilibrium();
        state.equilibrium().setControlVariablesP(p);
        state.equilibrium().setControlVariablesQ(q);
        state.equilibrium().setInputValues(w);
        state.equilibrium().setInitialComponentAmounts(c);

        auto const& pp = state.equilibrium().p();
        auto const& ww = state.equilibrium().w();

        const auto T = getT(pp, ww); // get temperature from predicted *p* or given *w*
        const auto P = getP(pp, ww); // get pressure from predicted *p* or given *w*

        state.setTemperature(T);
        state.setPressure(P);
    }

    /// Perform a first-order Taylor prediction of the chemical potential of a species at given conditions.
    auto speciesChemicalPotentialPredicted(Index i, VectorXdConstRef const& dw, VectorXdConstRef const& dc) const -> double
    {
        assert(i < Nn);

        const auto dudw0 = sensitivity0.dudw(); // The derivatives *du/dw* of the chemical properties of the chemical system wrt *w*.
        const auto dudc0 = sensitivity0.dudc(); // The derivatives *du/dc* of the chemical properties of the chemical system wrt *c*.

        const auto dmuidw0 = dudw0.row(Nu - Nn + i); // The derivatives *dμ[i]/dw* of the chemical potential of the i-th species.
        const auto dmuidc0 = dudc0.row(Nu - Nn + i); // The derivatives *dμ[i]/dc* of the chemical potential of the i-th species.
        const auto mui0 = u0[Nu - Nn + i];

        return mui0 + dmuidw0.dot(dw) + dmuidc0.dot(dc);
    }

    /// Return the chemical potential of a species at given reference conditions.
    auto speciesChemicalPotentialReference(Index i) const -> double
    {
        assert(i < Nn);
        return u0[Nu - Nn + i];
    }
};

EquilibriumPredictor::EquilibriumPredictor(ChemicalState const& state0, EquilibriumSensitivity const& sensitivity0)
: pimpl(new Impl(state0, sensitivity0))
{}

EquilibriumPredictor::EquilibriumPredictor(EquilibriumPredictor const& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumPredictor::~EquilibriumPredictor()
{}

auto EquilibriumPredictor::operator=(EquilibriumPredictor other) -> EquilibriumPredictor&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumPredictor::predict(ChemicalState& state, EquilibriumConditions const& conditions) const -> void
{
    pimpl->predict(state, conditions);
}

auto EquilibriumPredictor::speciesChemicalPotentialPredicted(Index ispecies, VectorXdConstRef const& dw, VectorXdConstRef const& dc) const -> double
{
    return pimpl->speciesChemicalPotentialPredicted(ispecies, dw, dc);
}

auto EquilibriumPredictor::speciesChemicalPotentialReference(Index ispecies) const -> double
{
    return pimpl->speciesChemicalPotentialReference(ispecies);
}

} // namespace Reaktoro
