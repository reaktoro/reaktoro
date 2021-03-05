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

#include "EquilibriumPredictor.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>

namespace Reaktoro {
namespace {

/// The auxiliary return type for the temperature and pressure getter functions below.
using GetterFn = Fn<double(VectorXdConstRef, VectorXdConstRef)>;

/// Create the lambda function that extracts temperature from either `p` or `w`.
auto getTemperatureFn(const Params& w) -> GetterFn
{
    const auto idxT = w.find("T");
    const auto knownT = idxT < w.size(); // T is known if in w (vector of input parameters)
    if(knownT) return [idxT](VectorXdConstRef p, VectorXdConstRef w) { return w[idxT]; };
    else return [](VectorXdConstRef p, VectorXdConstRef w) { return p[0]; };
}

/// Create the lambda function that extracts pressure from either `p` or `w`.
auto getPressureFn(const Params& w) -> GetterFn
{
    const auto idxP = w.find("P");
    const auto knownP = idxP < w.size(); // P is known if in w (vector of input parameters)
    if(knownP) return [idxP](VectorXdConstRef p, VectorXdConstRef w) { return w[idxP]; };
    const auto idxT = w.find("T");
    const auto knownT = idxT < w.size(); // T is known if in w (vector of input parameters)
    if(knownT) return [](VectorXdConstRef p, VectorXdConstRef w) { return p[0]; };
    else return [](VectorXdConstRef p, VectorXdConstRef w) { return p[1]; };
}

} // namespace

struct EquilibriumPredictor::Impl
{
    const VectorXd n0;    ///< The species amounts *n* at the reference equilibrium state.
    const VectorXd p0;    ///< The control variables *p* at the reference equilibrium state.
    const VectorXd w0;    ///< The input parameters *w* at the reference equilibrium state.
    const VectorXd b0;    ///< The component amounts *b* at the reference equilibrium state.
    const VectorXd u0;    ///< The chemical properties *u* at the reference equilibrium state.
    const MatrixXd dndw0; ///< The derivatives *dn/dw* at the reference equilibrium state
    const MatrixXd dpdw0; ///< The derivatives *dp/dw* at the reference equilibrium state
    const MatrixXd dudw0; ///< The derivatives *du/dw* at the reference equilibrium state
    const MatrixXd dndb0; ///< The derivatives *dn/db* at the reference equilibrium state
    const MatrixXd dpdb0; ///< The derivatives *dp/db* at the reference equilibrium state
    const MatrixXd dudb0; ///< The derivatives *du/db* at the reference equilibrium state
    VectorXd n;           ///< The species amounts *n* at the predicted equilibrium state.
    VectorXd p;           ///< The control variables *p* at the predicted equilibrium state.
    VectorXd w;           ///< The input parameters *w* at the predicted equilibrium state.
    VectorXd b;           ///< The component amounts *b* at the predicted equilibrium state.
    VectorXd u;           ///< The chemical properties *u* at the predicted equilibrium state.
    GetterFn getT;        ///< The function that gets temperature from either *p* or *w* depending if it is known or unwknon in the equilibrium calculation.
    GetterFn getP;        ///< The function that gets pressure from either *p* or *w* depending if it is known or unwknon in the equilibrium calculation.

    /// Construct a EquilibriumPredictor object.
    Impl(const ChemicalState& state0, const EquilibriumSensitivity& sensitivity0)
    : n0(state0.speciesAmounts()),
      p0(state0.equilibrium().p()),
      w0(state0.equilibrium().w()),
      b0(state0.equilibrium().b()),
      u0(state0.props()),
      dndw0(sensitivity0.dndw()),
      dpdw0(sensitivity0.dpdw()),
      dudw0(sensitivity0.dudw()),
      dndb0(sensitivity0.dndb()),
      dpdb0(sensitivity0.dpdb()),
      dudb0(sensitivity0.dudb()),
      getT(getTemperatureFn(state0.equilibrium().w())),
      getP(getPressureFn(state0.equilibrium().w()))
    {
    }

    /// Perform a first-order Taylor prediction of the chemical state at given conditions.
    auto predict(ChemicalState& state, const EquilibriumConditions& conditions) -> void
    {
        w = conditions.params();
        b = conditions.initialComponentAmounts();
        const auto dw = w - w0;
        const auto db = b - b0;
        n.noalias() = n0 + dndw0*dw + dndb0*db;
        p.noalias() = p0 + dpdw0*dw + dpdb0*db;
        u.noalias() = u0 + dudw0*dw + dudb0*db;
        const auto T = getT(p, w); // get temperature from predicted *p* or given *w*
        const auto P = getP(p, w); // get pressure from predicted *p* or given *w*
        state.setTemperature(T);
        state.setPressure(P);
        state.setSpeciesAmounts(n);
        state.props().update(u);
    }
};

EquilibriumPredictor::EquilibriumPredictor(const ChemicalState& state0, const EquilibriumSensitivity& sensitivity0)
: pimpl(new Impl(state0, sensitivity0))
{}

EquilibriumPredictor::EquilibriumPredictor(const EquilibriumPredictor& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumPredictor::~EquilibriumPredictor()
{}

auto EquilibriumPredictor::operator=(EquilibriumPredictor other) -> EquilibriumPredictor&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumPredictor::predict(ChemicalState& state, const EquilibriumConditions& conditions) -> void
{
    pimpl->predict(state, conditions);
}

} // namespace Reaktoro
