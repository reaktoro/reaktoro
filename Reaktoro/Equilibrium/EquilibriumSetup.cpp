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

/// Create the lambda function that extracts temperature value from either `p` or `params`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// temperature in unknown, then temperature should be extracted from vector `p`.
/// Otherwise, temperature is known and available in the Params object `params`.
/// @param specs The specifications of the equilibrium solver
auto createTemperatureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, const Params&)>
{
    const auto unknownT = specs.isTemperatureUnknown();

    if(unknownT)
        return [](VectorXrConstRef p, const Params& params) { return p[0]; };

    return [](VectorXrConstRef p, const Params& params) { return params.get("T").value(); };
}

/// Create the lambda function that extracts pressure value from either `p` or `params`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// pressure in unknown, then pressure should be extracted from vector `p`.
/// Otherwise, pressure is known and available in the Params object `params`.
/// @param specs The specifications of the equilibrium solver
auto createPressureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, const Params&)>
{
    const auto unknownT = specs.isTemperatureUnknown();
    const auto unknownP = specs.isPressureUnknown();

    if(unknownT && unknownP)
        return [](VectorXrConstRef p, const Params& params) { return p[1]; };

    if(unknownP)
        return [](VectorXrConstRef p, const Params& params) { return p[0]; };

    return [](VectorXrConstRef p, const Params& params) { return params.get("P").value(); };
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

    /// The temperature getter function for the given equilibrium specifications. @see createTemperatureGetterFn
    const Fn<real(VectorXrConstRef, const Params&)> getT;

    /// The pressure getter function for the given equilibrium specifications. @see createPressureGetterFn
    const Fn<real(VectorXrConstRef, const Params&)> getP;

    /// The auxiliary chemical properties of the system.
    ChemicalProps props;

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
    MatrixXd Hxc;  ///< The Jacobian of the objective gradient function with respect to params (Optima uses `c` for these parameter variables).
    VectorXr vp;   ///< The residuals of the equation constraints.
    MatrixXd Vpx;  ///< The Jacobian of the residuals of the equation constraints with respect to x = (n, q).
    MatrixXd Vpp;  ///< The Jacobian of the residuals of the equation constraints with respect to p.
    MatrixXd Vpc;  ///< The Jacobian of the residuals of the equation constraints with respect to params (Optima uses `c` for these parameter variables).
    Params params; ///< The current values of the input parameters for the equilibrium calculation (e.g, T, P, pH, pE, V, etc.).

    /// Construct an EquilibriumSetup::Impl object
    Impl(const EquilibriumSpecs& specs)
    : system(specs.system()), specs(specs), dims(specs), props(specs.system()),
      getT(createTemperatureGetterFn(specs)), getP(createPressureGetterFn(specs))
    {
        const auto Nn = system.species().size();
        const auto Np = numControlVariablesTypeP(specs);
        const auto Nq = numControlVariablesTypeQ(specs);
        const auto Nx = Nn + Nq;

        n.resize(Nn);
        p.resize(Np);
        q.resize(Nq);
        x.resize(Nx);
        gx.resize(Nx);
        Hxx.resize(Nx, Nx);
        Hxp.resize(Nx, Np);
        vp.resize(Np);
        Vpx.resize(Np, Nx);
        Vpp.resize(Np, Np);
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
        const auto Nc = numComponents(specs);
        const auto Nx = Nn + Nq;

        assert(Nc == Ne); // TODO: Remove this when EquilibriumReactions is implemented and inert reactions can be set

        MatrixXd Aex = zeros(Nc, Nx);

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
        const auto Nc = numComponents(specs);

        assert(Nc == Ne); // TODO: Remove this when EquilibriumReactions is implemented and inert reactions can be set

        MatrixXd Aep = zeros(Nc, Np);

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
        const auto Wn = system.formulaMatrix();
        const auto n0 = state0.speciesAmounts();
        const auto be = Wn * n0.matrix();
        return be;
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
    auto updateChemicalProps(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> void
    {
        //======================================================================
        // Update temperature (T) and pressure (P)
        //----------------------------------------------------------------------
        // If temperature or pressure are unknowns, get their current values
        // from vector p. If not, they have been given in the params object.
        //======================================================================
        T = getT(p, params);
        P = getP(p, params);

        assert(T > 0.0); // check if a proper value for T was given in p or params
        assert(P > 0.0); // check if a proper value for P was given in p or params

        //======================================================================
        // Update chemical properties of the chemical system
        //======================================================================
        const auto Nn = system.species().size();
        const auto n = x.head(Nn);

        props.update(T, P, n);
    }

    auto evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> real
    {
        updateChemicalProps(x, p, params);
        const auto& T = props.temperature();
        const auto& n = props.speciesAmounts();
        const auto& u = props.chemicalPotentials();
        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;
        const auto barrier = -tau * n.log().sum();
        return (n * u).sum()/RT + barrier; // the current Gibbs energy of the system (normalized by RT)
    }

    auto evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef
    {
        updateChemicalProps(x, p, params);

        const auto Nn = system.species().size();
        const auto Nq = numControlVariablesTypeQ(specs);

        const auto& T = props.temperature();
        const auto& n = props.speciesAmounts();
        const auto& u = props.chemicalPotentials();

        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;

        auto gn = gx.head(Nn); // where we set the chemical potentials of the species
        auto gq = gx.tail(Nq); // where we set the desired chemical potentials of some substances

        gn = u/RT - tau/n; // set the current chemical potentials of species (normalized by RT)

        const auto& uconstraints = specs.constraintsChemicalPotentialType();

        for(auto i = 0; i < Nq; ++i)
            gq[i] = uconstraints[i].fn(props)/RT;

        return gx;
    }

    auto evalObjectiveHessianX(VectorXrConstRef xconst, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
    {
        x = xconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalObjectiveGradX(x, p, params);
        };
        Hxx = jacobian(fn, wrt(x), at(x, p, params));
        return Hxx;
    }

    auto evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef pconst, const Params& params) -> MatrixXdConstRef
    {
        p = pconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalObjectiveGradX(x, p, params);
        };
        Hxp = jacobian(fn, wrt(p), at(x, p, params));
        return Hxp;
    }

    auto evalObjectiveHessianParams(VectorXrConstRef x, VectorXrConstRef p, const Params& paramsconst) -> MatrixXdConstRef
    {
        params = paramsconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalObjectiveGradX(x, p, params);
        };
        Hxc = jacobian(fn, wrt(params), at(x, p, params));
        return Hxc;
    }

    auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef
    {
        updateChemicalProps(x, p, params);

        const auto Np = numControlVariablesTypeP(specs);

        const auto& econstraints = specs.constraintsEquationType();

        for(auto i = 0; i < Np; ++i)
            vp[i] = econstraints[i].fn(props);

        return vp;
    }

    auto evalEquationConstraintsGradX(VectorXrConstRef xconst, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
    {
        x = xconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalEquationConstraints(x, p, params);
        };
        Vpx = jacobian(fn, wrt(x), at(x, p, params));
        return Vpx;
    }

    auto evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef pconst, const Params& params) -> MatrixXdConstRef
    {
        p = pconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalEquationConstraints(x, p, params);
        };
        Vpp = jacobian(fn, wrt(p), at(x, p, params));
        return Vpp;
    }

    auto evalEquationConstraintsGradParams(VectorXrConstRef x, VectorXrConstRef p, const Params& paramsconst) -> MatrixXdConstRef
    {
        params = paramsconst;
        auto fn = [&](VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXr
        {
            return evalEquationConstraints(x, p, params);
        };
        Vpc = jacobian(fn, wrt(params), at(x, p, params));
        return Vpc;
    }

    auto extractTemperature(VectorXrConstRef p, const Params& params) const -> real
    {
        return getT(p, params);
    }

    auto extractPressure(VectorXrConstRef p, const Params& params) const -> real
    {
        return getP(p, params);
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

auto EquilibriumSetup::evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> real
{
    return pimpl->evalObjectiveValue(x, p, params);
}

auto EquilibriumSetup::evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef
{
    return pimpl->evalObjectiveGradX(x, p, params);
}

auto EquilibriumSetup::evalObjectiveHessianX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianX(x, p, params);
}

auto EquilibriumSetup::evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianP(x, p, params);
}

auto EquilibriumSetup::evalObjectiveHessianParams(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalObjectiveHessianParams(x, p, params);
}

auto EquilibriumSetup::evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef
{
    return pimpl->evalEquationConstraints(x, p, params);
}

auto EquilibriumSetup::evalEquationConstraintsGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradX(x, p, params);
}

auto EquilibriumSetup::evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradP(x, p, params);
}

auto EquilibriumSetup::evalEquationConstraintsGradParams(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef
{
    return pimpl->evalEquationConstraintsGradParams(x, p, params);
}

auto EquilibriumSetup::extractTemperature(VectorXrConstRef p, const Params& params) const -> real
{
    return pimpl->extractTemperature(p, params);
}

auto EquilibriumSetup::extractPressure(VectorXrConstRef p, const Params& params) const -> real
{
    return pimpl->extractPressure(p, params);
}

} // namespace Reaktoro
