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

#include "EquilibriumSetup.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {
namespace {

/// Return the number of control variables *p* in the equilibrium specifications.
/// The control variables *p* are temperature, pressure and/or amounts of explicit titrants.
auto numControlVariablesTypeP(const EquilibriumSpecs& specs) -> Index
{
    return specs.numControlVariables() - specs.numTitrantsImplicit();
}

/// Return the number of control variables *q* in the equilibrium specifications.
/// The control variables *q* are the amounts of implicit titrants.
auto numControlVariablesTypeQ(const EquilibriumSpecs& specs) -> Index
{
    return specs.numTitrantsImplicit();
}

/// Return the number of components in the equilibrium specifications.
/// These are the independent primary species in the equilibrium problem.
auto numComponents(const EquilibriumSpecs& specs) -> Index
{
    return specs.system().elements().size() + 1;
}

} // namespace

using autodiff::jacobian;
using autodiff::wrt;
using autodiff::at;

struct EquilibriumSetup::Impl
{
    /// The chemical system associated with this equilibrium problem
    const ChemicalSystem system;

    /// The equilibrium specifications associated with this equilibrium problem
    const EquilibriumSpecs specs;

    /// The dimensions of the variables in the equilibrium problem.
    const EquilibriumDims dims;

    /// The auxiliary chemical properties of the system.
    EquilibriumProps props;

    /// The options for the solution of the equilibrium problem.
    EquilibriumOptions options;

    real T = 0.0;  ///< The current temperature of the chemical system (in K).
    real P = 0.0;  ///< The current pressure of the chemical system (in Pa).
    ArrayXr n;     ///< The current amounts of the species (in mol).
    ArrayXr p;     ///< The current values of the introduced *p* control variables (temperature, pressure, amounts of explicit titrants).
    ArrayXr q;     ///< The current values of the introduced *q* control variables (amounts of implicit titrants whose chemical potentials are constrained).
    VectorXr x;    ///< The auxiliary vector x = (n, q).
    VectorXr gx;   ///< The gradient of the objective function with respect to x = (n, q).
    MatrixXd Hxx;  ///< The Jacobian of the objective gradient function with respect to x = (n, q).
    MatrixXd Hxp;  ///< The Jacobian of the objective gradient function with respect to p.
    MatrixXd Hxc;  ///< The Jacobian of the objective gradient function with respect to w (Optima uses `c` for these parameter variables).
    VectorXr vp;   ///< The residuals of the equation constraints.
    MatrixXd Vpx;  ///< The Jacobian of the residuals of the equation constraints with respect to x = (n, q).
    MatrixXd Vpp;  ///< The Jacobian of the residuals of the equation constraints with respect to p.
    MatrixXd Vpc;  ///< The Jacobian of the residuals of the equation constraints with respect to w (Optima uses `c` for these parameter variables).
    Params w; ///< The current values of the input parameters for the equilibrium calculation (e.g, T, P, pH, pE, V, etc.).

    /// Construct an EquilibriumSetup::Impl object
    Impl(const EquilibriumSpecs& specs)
    : system(specs.system()), specs(specs), dims(specs), props(specs),
      w(specs.params())
    {
        const auto Nn = system.species().size();
        const auto Np = numControlVariablesTypeP(specs);
        const auto Nq = numControlVariablesTypeQ(specs);
        const auto Nx = Nn + Nq;
        const auto Nb = numComponents(specs);
        const auto Nc = w.size() + Nb;

        n.resize(Nn);
        p.resize(Np);
        q.resize(Nq);
        x.resize(Nx);
        gx.resize(Nx);
        Hxx.resize(Nx, Nx);
        Hxp.resize(Nx, Np);
        Hxc.resize(Nx, Nc);
        vp.resize(Np);
        Vpx.resize(Np, Nx);
        Vpp.resize(Np, Np);
        Vpc.resize(Np, Nc);
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) const -> void
    {
        const auto Ne = system.elements().size() + 1;
        assert(vec.size() == Ne);
        vec.fill(0.0);
        vec[Ne - 1] = formula.charge(); // last entry in the column vector is charge of substance
        for(const auto& [element, coeff] : formula.elements()) {
            const auto ielem = system.elements().index(element);
            vec[ielem] = coeff;
        }
    }

    /// Assemble the coefficient matrix `Aex` in optimization problem.
    auto assembleMatrixAex() const -> MatrixXd
    {
        const auto Ne = system.elements().size() + 1;
        const auto Nn = system.species().size();
        const auto Nq = numControlVariablesTypeQ(specs);
        const auto Nb = numComponents(specs);
        const auto Nx = Nn + Nq;

        assert(Nb == Ne); // TODO: Remove this when EquilibriumReactions is implemented and inert reactions can be set

        MatrixXd Aex = zeros(Nb, Nx);

        auto Wn = Aex.topLeftCorner(Ne, Nn);  // the formula matrix of the species
        auto Wq = Aex.topRightCorner(Ne, Nq); // the formula matrix of the implicit titrants

        Wn = system.formulaMatrix();

        auto j = 0;
        for(const auto& formula : specs.titrantsImplicit())
            assembleFormulaVector(Wq.col(j++), formula);

        return Aex;
    }

    /// Assemble the coefficient matrix `Aep` in optimization problem.
    auto assembleMatrixAep() const -> MatrixXd
    {
        const auto Ne = system.elements().size() + 1;
        const auto Np = numControlVariablesTypeP(specs);
        const auto Nb = numComponents(specs);

        assert(Nb == Ne); // TODO: Remove this when EquilibriumReactions is implemented and inert reactions can be set

        MatrixXd Aep = zeros(Nb, Np);

        auto Wp = Aep.topRows(Ne); // the formula matrix of temperature, pressure and explicit titrants

        auto j = specs.isTemperatureUnknown() + specs.isPressureUnknown(); // skip columns corresponding to T and P in p (if applicable), since these are zeros
        for(const auto& formula : specs.titrantsExplicit())
            assembleFormulaVector(Wp.col(j++), formula);

        return Aep;
    }

    /// Assemble the right-hand side vector `be` in the optimization problem.
    auto assembleVectorBe(const EquilibriumConditions& conditions, const ChemicalState& state0) -> VectorXr
    {
        if(conditions.initialComponentAmounts().size())
            return conditions.initialComponentAmounts();
        else return conditions.initialComponentAmountsCompute(state0.speciesAmounts());
    }

    /// Assemble the lower bound vector `xlower` in the optimization problem where *x = (n, q)*.
    auto assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
    {
        const auto Nn = system.species().size();
        const auto Nq = numControlVariablesTypeQ(specs);
        const auto Nx = Nn + Nq;
        VectorXd xlower = constants(Nx, -inf);
        auto nlower = xlower.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.speciesCannotDecreaseBelow()) nlower[i] = val;
        for(auto i : restrictions.speciesCannotDecrease()) nlower[i] = n0[i]; // this comes after, in case a species cannot strictly decrease
        for(auto& val : nlower) val = std::max(val, options.epsilon); // ensure the upper bounds of the species amounts are not below the minimum amount value given in EquilibriumOptions::epsilon. TODO: Issue a warning when lower/upper bound of a species amount is changed to EquilibriumOptions::epsilon.
        return xlower;
    }

    /// Assemble the upper bound vector `xupper` in the optimization problem where *x = (n, q)*.
    auto assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
    {
        const auto Nn = system.species().size();
        const auto Nq = numControlVariablesTypeQ(specs);
        const auto Nx = Nn + Nq;
        VectorXd xupper = constants(Nx, inf);
        auto nupper = xupper.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.speciesCannotIncreaseAbove()) nupper[i] = val;
        for(auto i : restrictions.speciesCannotIncrease()) nupper[i] = n0[i]; // this comes after, in case a species cannot strictly increase
        for(auto& val : nupper) val = std::max(val, options.epsilon); // ensure the upper bounds of the species amounts are not below the minimum amount value given in EquilibriumOptions::epsilon.
        return xupper;
    }

    /// Update the chemical properties of the system.
    auto updateChemicalProps(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> void
    {
        const auto n = x.head(dims.Nn);
        props.update(n, p, w);
    }

    auto evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> real
    {
        updateChemicalProps(x, p, w);
        const auto& cprops = props.chemicalProps();
        const auto& T = cprops.temperature();
        const auto& n = cprops.speciesAmounts();
        const auto& u = cprops.chemicalPotentials();
        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;
        const auto barrier = -tau * n.log().sum();
        return (n * u).sum()/RT + barrier; // the current Gibbs energy of the system (normalized by RT)
    }

    auto evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXrConstRef
    {
        updateChemicalProps(x, p, w);

        const auto Nn = system.species().size();
        const auto Nq = numControlVariablesTypeQ(specs);

        const auto& cprops = props.chemicalProps();
        const auto& T = cprops.temperature();
        const auto& n = cprops.speciesAmounts();
        const auto& u = cprops.chemicalPotentials();

        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;

        auto gn = gx.head(Nn); // where we set the chemical potentials of the species
        auto gq = gx.tail(Nq); // where we set the desired chemical potentials of some substances

        gn = u/RT - tau/n; // set the current chemical potentials of species (normalized by RT)

        const auto& uconstraints = specs.constraintsChemicalPotentialType();

        for(auto i = 0; i < Nq; ++i)
            gq[i] = uconstraints[i].fn(cprops)/RT;

        return gx;
    }

    auto evalObjectiveHessianX(VectorXrConstRef xconst, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
    {
        x = xconst;
        auto fn = [this](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalObjectiveGradX(x, p, w);
        };
        Hxx = jacobian(fn, wrt(x), at(x, p, w));
        return Hxx;
    }

    auto evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef pconst, const Params& w) -> MatrixXdConstRef
    {
        p = pconst;
        auto fn = [this](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalObjectiveGradX(x, p, w);
        };
        Hxp = jacobian(fn, wrt(p), at(x, p, w));
        return Hxp;
    }

    auto evalObjectiveHessianParams(VectorXrConstRef x, VectorXrConstRef p, const Params& wconst) -> MatrixXdConstRef
    {
        w = wconst;
        auto fn = [this](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalObjectiveGradX(x, p, w);
        };
        Hxc.leftCols(w.size()) = jacobian(fn, wrt(w), at(x, p, w));
        Hxc.rightCols(dims.Nb).setZero(); // these are derivatives wrt amounts of components
        return Hxc;
    }

    auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXrConstRef
    {
        updateChemicalProps(x, p, w);

        const auto Np = numControlVariablesTypeP(specs);

        const auto& econstraints = specs.constraintsEquationType();
        const auto& cprops = props.chemicalProps();

        for(auto i = 0; i < Np; ++i)
            vp[i] = econstraints[i].fn(cprops);

        return vp;
    }

    auto evalEquationConstraintsGradX(VectorXrConstRef xconst, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
    {
        x = xconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalEquationConstraints(x, p, w);
        };
        Vpx = jacobian(fn, wrt(x), at(x, p, w));
        return Vpx;
    }

    auto evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef pconst, const Params& w) -> MatrixXdConstRef
    {
        p = pconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalEquationConstraints(x, p, w);
        };
        Vpp = jacobian(fn, wrt(p), at(x, p, w));
        return Vpp;
    }

    auto evalEquationConstraintsGradParams(VectorXrConstRef x, VectorXrConstRef p, const Params& wconst) -> MatrixXdConstRef
    {
        w = wconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            return evalEquationConstraints(x, p, w);
        };
        Vpc.leftCols(dims.Nw) = jacobian(fn, wrt(w), at(x, p, w));
        Vpc.rightCols(dims.Nb).setZero(); // these are derivatives wrt amounts of components, which are zero
        return Vpc;
    }
};

EquilibriumSetup::EquilibriumSetup(const EquilibriumSpecs& conditions)
: pimpl(new Impl(conditions))
{}

EquilibriumSetup::EquilibriumSetup(const EquilibriumSetup& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSetup::~EquilibriumSetup()
{}

auto EquilibriumSetup::operator=(EquilibriumSetup other) -> EquilibriumSetup&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSetup::setOptions(const EquilibriumOptions& opts) -> void
{
    pimpl->options = opts;
}

auto EquilibriumSetup::dims() const -> const EquilibriumDims&
{
    return pimpl->dims;
}

auto EquilibriumSetup::options() const -> const EquilibriumOptions&
{
    return pimpl->options;
}

auto EquilibriumSetup::assembleMatrixAex() const -> MatrixXd
{
    return pimpl->assembleMatrixAex();
}

auto EquilibriumSetup::assembleMatrixAep() const -> MatrixXd
{
    return pimpl->assembleMatrixAep();
}

auto EquilibriumSetup::assembleVectorBe(const EquilibriumConditions& conditions, const ChemicalState& state0) -> VectorXr
{
    return pimpl->assembleVectorBe(conditions, state0);
}

auto EquilibriumSetup::assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
{
    return pimpl->assembleLowerBoundsVector(restrictions, state0);
}

auto EquilibriumSetup::assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
{
    return pimpl->assembleUpperBoundsVector(restrictions, state0);
}

auto EquilibriumSetup::evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> real
{
    return pimpl->evalObjectiveValue(x, p, w);
}

auto EquilibriumSetup::evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXrConstRef
{
    return pimpl->evalObjectiveGradX(x, p, w);
}

auto EquilibriumSetup::evalObjectiveHessianX(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianX(x, p, w);
}

auto EquilibriumSetup::evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianP(x, p, w);
}

auto EquilibriumSetup::evalObjectiveHessianParams(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianParams(x, p, w);
}

auto EquilibriumSetup::evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> VectorXrConstRef
{
    return pimpl->evalEquationConstraints(x, p, w);
}

auto EquilibriumSetup::evalEquationConstraintsGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradX(x, p, w);
}

auto EquilibriumSetup::evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradP(x, p, w);
}

auto EquilibriumSetup::evalEquationConstraintsGradParams(VectorXrConstRef x, VectorXrConstRef p, const Params& w) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradParams(x, p, w);
}

auto EquilibriumSetup::assembleChemicalPropsJacobianBegin() -> void
{
    pimpl->props.assembleFullJacobianBegin();
}

auto EquilibriumSetup::assembleChemicalPropsJacobianEnd() -> void
{
    pimpl->props.assembleFullJacobianEnd();
}

auto EquilibriumSetup::equilibriumProps() const -> const EquilibriumProps&
{
    return pimpl->props;
}

auto EquilibriumSetup::chemicalProps() const -> const ChemicalProps&
{
    return pimpl->props.chemicalProps();
}

} // namespace Reaktoro
