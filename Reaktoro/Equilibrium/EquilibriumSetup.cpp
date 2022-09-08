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

#include "EquilibriumSetup.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {
namespace {

/// Assemble the vector with the element and charge coefficients of a chemical formula.
auto assembleFormulaVector(ChemicalFormula const& formula, ChemicalSystem const& system) -> VectorXd
{
    const auto Ne = system.elements().size();
    VectorXd res = zeros(Ne + 1);
    res[Ne] = formula.charge(); // last entry in the column vector is charge of substance
    for(auto const& [element, coeff] : formula.elements()) {
        const auto ielem = system.elements().index(element);
        res[ielem] = coeff;
    }
    return res;
}

/// Assemble the coefficient matrix `Aex` in optimization problem.
auto assembleMatrixAex(EquilibriumSpecs const& specs) -> MatrixXd
{
    const EquilibriumDims dims(specs);

    MatrixXd Aex = zeros(dims.Nb, dims.Nx);

    auto Wn = Aex.topLeftCorner(dims.Ne + 1, dims.Nn);  // the formula matrix of the species with respect to elements and charge
    auto Wq = Aex.topRightCorner(dims.Ne + 1, dims.Nq); // the formula matrix of the implicit titrants with respect to elements and charge
    auto Kn = Aex.bottomLeftCorner(dims.Nr, dims.Nn);   // the coefficient matrix of the reactivity constraints with respect to the species amount variables

    Wn = specs.system().formulaMatrix();

    for(auto [i, formula] : enumerate(specs.titrantsImplicit()))
        Wq.col(i) = assembleFormulaVector(formula, specs.system());

    for(auto [i, rconstraint] : enumerate(specs.reactivityConstraints()))
        Kn.row(i) = rconstraint.Kn;

    return Aex;
}

/// Assemble the coefficient matrix `Aep` in optimization problem.
auto assembleMatrixAep(EquilibriumSpecs const& specs) -> MatrixXd
{
    const EquilibriumDims dims(specs);

    MatrixXd Aep = zeros(dims.Nb, dims.Np);

    auto Wp = Aep.topRows(dims.Ne + 1); // the formula matrix of the p variables with respect to elements and charge (e.g., temperature, pressure, custom variables, and explicit titrants)
    auto Kp = Aep.bottomRows(dims.Nr);  // the coefficient matrix of the reactivity constraints with respect to the p control variables

    auto offset = specs.isTemperatureUnknown() + specs.isPressureUnknown(); // skip columns corresponding to T and P variables in p (if applicable, if they are unknown), since these columns are zeros

    for(auto [i, formula] : enumerate(specs.titrantsExplicit()))
        Wp.col(i + offset) = assembleFormulaVector(formula, specs.system());

    for(auto [i, rconstraint] : enumerate(specs.reactivityConstraints()))
        if(rconstraint.Kp.size())
            Kp.row(i) = rconstraint.Kp;

    return Aep;
}

} // namespace

struct EquilibriumSetup::Impl
{
    const ChemicalSystem system;  ///< The chemical system associated with this equilibrium problem.
    const EquilibriumSpecs specs; ///< The equilibrium specifications associated with this equilibrium problem.
    const EquilibriumDims dims;   ///< The dimensions of the variables in the equilibrium problem.
    const MatrixXd Aex;           ///< The coefficient matrix Aex in the optimization problem.
    const MatrixXd Aep;           ///< The coefficient matrix Aep in the optimization problem.
    EquilibriumProps props;       ///< The auxiliary chemical properties of the system.
    EquilibriumOptions options;   ///< The options for the solution of the equilibrium problem.
    VectorXr x;                   ///< The auxiliary vector x = (n, q).
    VectorXr n;                   ///< The current amounts of the species (in mol).
    VectorXr q;                   ///< The current values of the introduced *q* control variables (amounts of implicit titrants whose chemical potentials are constrained).
    VectorXr p;                   ///< The current values of the introduced *p* control variables (temperature, pressure, amounts of explicit titrants).
    VectorXr w;                   ///< The current values of the input variables for the equilibrium calculation (e.g, T, P, pH, pE, V, etc.).
    VectorXr F;                   ///< The auxiliary vector F = (gx, vp) containing first-order optimality conditions (gx) and residuals of equation constraints (vp)
    real f;                       ///< The objective function value (Gibbs energy of the system).
    VectorXd gx;                  ///< The gradient of the objective function with respect to x = (n, q).
    VectorXd vp;                  ///< The residual vector of the equilibrium equation constraints (of size p).
    MatrixXd Hxx;                 ///< The Jacobian of the objective gradient function with respect to x = (n, q).
    MatrixXd Hxp;                 ///< The Jacobian of the objective gradient function with respect to p.
    MatrixXd Hxc;                 ///< The Jacobian of the objective gradient function with respect to c = (w, b).
    MatrixXd Vpx;                 ///< The Jacobian of vp with respect to x = (n, q).
    MatrixXd Vpp;                 ///< The Jacobian of vp with respect to p.
    MatrixXd Vpc;                 ///< The Jacobian of vp with respect to c = (w, b).
    ArrayXr mu;                   ///< The auxiliary vector of chemical potentials of the species.
    VectorXl isbasicvar;          /// The bitmap that indicates which variables in x = (n, q) are currently basic variables.
    Indices ipps;                 ///< The indices of the pure phase species (i.e., species composing single-phase species, whose chemical potentials do not depend on composition)

    // -------------------------------------------- //
    // ------ CONVENIENT AUXILIARY VARIABLES ------ //
    // -------------------------------------------- //

    const Index Nb; ///< The number of conservative components
    const Index Ne; ///< The number of elements
    const Index Nn; ///< The number of chemical species
    const Index Np; ///< The number of p control variables
    const Index Nq; ///< The number of q control variables
    const Index Nx; ///< The number of x variables with x = (n, q)
    const Index Nw; ///< The number of input variables w in the chemical equilibrium problem
    const Index Nc; ///< The number of input variables in c = (w, b)
    const Index Nr; ///< The number of reactivity constraints (aka restricted reactions)

    /// Construct an EquilibriumSetup::Impl object
    Impl(EquilibriumSpecs const& specs)
    : system(specs.system()), specs(specs), dims(specs), Aex(assembleMatrixAex(specs)), Aep(assembleMatrixAep(specs)), props(specs),
      Nb(dims.Nb), Ne(dims.Ne), Nn(dims.Nn), Np(dims.Np), Nq(dims.Nq), Nx(dims.Nx), Nw(dims.Nw), Nc(dims.Nw + dims.Nb), Nr(dims.Nr)
    {
        x.resize(Nx);
        n.resize(Nn);
        q.resize(Nq);
        p.resize(Np);
        w.resize(Nw);
        F.resize(Nx + Np);
        gx.resize(Nx);
        vp.resize(Np);
        Hxx.resize(Nx, Nx);
        Hxp.resize(Nx, Np);
        Hxc.resize(Nx, Nc);
        Vpx.resize(Np, Nx);
        Vpp.resize(Np, Np);
        Vpc.resize(Np, Nc);
        mu.resize(Nn);

        isbasicvar.resize(Nx);

        // Initialize the indices of the pure phase species
        auto offset = 0;
        for(auto const& phase : system.phases())
        {
            const auto size = phase.species().size();
            if(size == 1)
                ipps.push_back(offset);
            offset += size;
        }
    }

    auto assembleLowerBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd
    {
        VectorXd xlower = constants(Nx, -inf);
        auto nlower = xlower.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.speciesCannotDecreaseBelow()) nlower[i] = val;
        for(auto i : restrictions.speciesCannotDecrease()) nlower[i] = n0[i]; // this comes after, in case a species cannot strictly decrease
        for(auto& val : nlower) val = std::max(val, options.epsilon); // ensure the upper bounds of the species amounts are not below the minimum amount value given in EquilibriumOptions::epsilon. TODO: Issue a warning when lower/upper bound of a species amount is changed to EquilibriumOptions::epsilon.
        return xlower;
    }

    auto assembleUpperBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd
    {
        VectorXd xupper = constants(Nx, inf);
        auto nupper = xupper.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.speciesCannotIncreaseAbove()) nupper[i] = val;
        for(auto i : restrictions.speciesCannotIncrease()) nupper[i] = n0[i]; // this comes after, in case a species cannot strictly increase
        for(auto& val : nupper) val = std::max(val, options.epsilon); // ensure the upper bounds of the species amounts are not below the minimum amount value given in EquilibriumOptions::epsilon.
        return xupper;
    }

    auto update(VectorXrConstRef xx, VectorXrConstRef pp, VectorXrConstRef ww) -> void
    {
        x = xx;
        n = xx.head(Nn);
        q = xx.tail(Nq);
        p = pp;
        w = ww;

        props.update(n, p, w, options.use_ideal_activity_models);

        updateF();
        updateGibbsEnergy(); // let this after updateF because of update in mu performed by updateF
        gx = F.head(Nx);
        vp = F.tail(Np);
    }

    auto updateGradX(VectorXlConstRef ibasicvars) -> void
    {
        isbasicvar.fill(false);
        isbasicvar(ibasicvars).fill(true);

        // Update Hxx and Vpx
        for(auto i = 0; i < Nn; ++i)
        {
            updateFx(i);
            Hxx.col(i) = grad(F.head(Nx));
            Vpx.col(i) = grad(F.tail(Np));
        }
        Hxx.rightCols(Nq).fill(0.0); // these are derivatives w.r.t. amounts of implicit titrants q
        Vpx.rightCols(Nq).fill(0.0); // these are derivatives w.r.t. amounts of implicit titrants q
    }

    auto updateGradP() -> void
    {
        // Update Hxp and Vpp
        for(auto i = 0; i < Np; ++i)
        {
            updateFp(i);
            Hxp.col(i) = grad(F.head(Nx));
            Vpp.col(i) = grad(F.tail(Np));
        }
    }

    auto updateGradW() -> void
    {
        // Update Hxc and Vpc
        for(auto i = 0; i < Nw; ++i)
        {
            updateFw(i);
            Hxc.col(i) = grad(F.head(Nx));
            Vpc.col(i) = grad(F.tail(Np));
        }
        Hxc.rightCols(Nb).fill(0.0); // these are derivatives w.r.t. amounts of conservative components
        Vpc.rightCols(Nb).fill(0.0); // these are derivatives w.r.t. amounts of conservative components
    }

    auto updateF() -> void
    {
        auto const& qvars = specs.controlVariablesQ();
        auto const& pvars = specs.controlVariablesP();

        auto const& econstraints = specs.equationConstraints();

        auto const& state = props.chemicalState();

        auto const& T = state.temperature();
        auto const& n = state.speciesAmounts();

        // Update the vector of species chemical potentials in case there are p variables associated to them
        mu = state.props().speciesChemicalPotentials();
        for(auto i = 0; i < Np; ++i)
            if(pvars[i].ispecies != Index(-1))
                mu[pvars[i].ispecies] = pvars[i].fn(state, p[i]);

        const auto RT  = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;

        auto gn = F.head(Nn);        // the segment in F where we set the chemical potentials of the species
        auto gq = F.segment(Nn, Nq); // the segment in F where we set the desired chemical potentials of some substances
        auto vp = F.tail(Np);        // the segment in F where we set the residuals of the equation constraints

        gn = mu/RT; // set the current chemical potentials of species (normalized by RT)
        gn(ipps).array() -= tau/n(ipps); // add log barrier contribution to pure phase species

        for(auto i = 0; i < Nq; ++i)
            gq[i] = qvars[i].fn(state, p, w)/RT;

        for(auto i = 0; i < Np; ++i)
            vp[i] = econstraints[i].fn(state, p, w);
    }

    auto updateFn(Index i) -> void
    {
        const auto useIdealModel = useIdealModelForGradWrtVariableN(i); // in case of little or no dependency of the thermochemical properties on n[i] (i.e., chemical props should have very little dependency in general on tiny species amounts)
        const auto inpw = i; // the index of n[i] in the extended vector (n, p, w)
        autodiff::seed(n[i]);
        props.update(n, p, w, useIdealModel, inpw);
        updateF();
        autodiff::unseed(n[i]);
    }

    auto updateFq(Index i) -> void
    {
        const auto useIdealModel = useIdealModelForGradWrtVariableQ(i); // in case of little or no dependency of the thermochemical properties on q[i] (i.e., chemical props has no dependency on amounts of implicit titrants such as [H+] when fixing pH)
        const auto inpw = -1; // the index of q[i] in the extended vector (n, p, w) is not defined
        autodiff::seed(q[i]);
        props.update(n, p, w, useIdealModel, inpw);
        updateF();
        autodiff::unseed(q[i]);
    }

    auto updateFp(Index i) -> void
    {
        assert(i < Np);
        const auto useIdealModel = useIdealModelForGradWrtVariableP(i); // in case of little or no dependency of the thermochemical properties on p[i] (e.g., chemical props has no dependency on the amount of a titrant, but it has on temperature and pressure if one of these are unknown p variables)
        const auto inpw = Nn + i; // the index of p[i] in the extended vector (n, p, w)
        autodiff::seed(p[i]);
        props.update(n, p, w, useIdealModel, inpw);
        updateF();
        autodiff::unseed(p[i]);
    }

    auto updateFw(Index i) -> void
    {
        assert(i < Nw);
        const auto useIdealModel = useIdealModelForGradWrtVariableW(i); // in case of little or no dependency of the thermochemical properties on w[i] (e.g., chemical props has no dependency on the designated value of pH, but it has on given values of temperature and pressure)
        const auto inpw = Nn + Np + i; // the index of w[i] in the extended vector (n, p, w)
        autodiff::seed(w[i]);
        props.update(n, p, w, useIdealModel, inpw);
        updateF();
        autodiff::unseed(w[i]);
    }

    auto updateFx(Index i) -> void
    {
        assert(i < Nx);
        if(i < Nn) updateFn(i);
        else updateFq(i);
    }

    auto updateGibbsEnergy() -> void
    {
        auto const& state = props.chemicalState();
        auto const& T = state.temperature();
        auto const& n = state.speciesAmounts();
        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;
        const auto barrier = -tau * n(ipps).log().sum();
        f = (n * mu).sum()/RT + barrier; // the current Gibbs energy of the system (normalized by RT)
    }

    auto getGibbsEnergy() -> real
    {
        return f;
    }

    auto getGibbsGradX() -> VectorXdConstRef
    {
        return gx;
    }

    auto getGibbsHessianX() -> MatrixXdConstRef
    {
        return Hxx;
    }

    auto getGibbsHessianP() -> MatrixXdConstRef
    {
        return Hxp;
    }

    auto getGibbsHessianC() -> MatrixXdConstRef
    {
        return Hxc;
    }

    auto getConstraintResiduals() -> VectorXdConstRef
    {
        return vp;
    }

    auto getConstraintResidualsGradX() -> MatrixXdConstRef
    {
        return Vpx;
    }

    auto getConstraintResidualsGradP() -> MatrixXdConstRef
    {
        return Vpp;
    }

    auto getConstraintResidualsGradC() -> MatrixXdConstRef
    {
        return Vpc;
    }

    auto usingPartiallyExactDerivatives() -> bool
    {
        return options.hessian == GibbsHessian::PartiallyExact;
    }

    auto usingDiagonalApproxDerivatives() -> bool
    {
        return options.hessian == GibbsHessian::ApproxDiagonal;
    }

    auto useIdealModelForGradWrtVariableN(Index i) -> bool
    {
        if(options.use_ideal_activity_models)
            return true;

        switch(options.hessian)
        {
        case GibbsHessian::Exact:          return false;
        case GibbsHessian::Approx:         return true;
        case GibbsHessian::ApproxDiagonal: return true;
        case GibbsHessian::PartiallyExact: return !isbasicvar[i];
        default:                           return !isbasicvar[i];
        }
    }

    auto useIdealModelForGradWrtVariableQ(Index i) -> bool
    {
        // Note: This should always be true because chemical properties do not
        // depend on the q variables (amounts of implicit titrants).
        return true;
    }

    auto useIdealModelForGradWrtVariableP(Index i) -> bool
    {
        // TODO: Improve efficiency by using ideal thermodynamic models when
        // computing derivatives with respect to some variables in p that play
        // no role in chemical properties of the system. For example, amounts
        // of explicit titrants.
        return false;
    }

    auto useIdealModelForGradWrtVariableW(Index i) -> bool
    {
        // TODO: Improve efficiency by using ideal thermodynamic models when
        // computing derivatives with respect to some variables in w that play
        // no role in chemical properties of the system. For example, pH input
        // value.
        return false;
    }
};

EquilibriumSetup::EquilibriumSetup(EquilibriumSpecs const& conditions)
: pimpl(new Impl(conditions))
{}

EquilibriumSetup::EquilibriumSetup(EquilibriumSetup const& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSetup::~EquilibriumSetup()
{}

auto EquilibriumSetup::operator=(EquilibriumSetup other) -> EquilibriumSetup&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSetup::setOptions(EquilibriumOptions const& opts) -> void
{
    pimpl->options = opts;
}

auto EquilibriumSetup::dims() const -> EquilibriumDims const&
{
    return pimpl->dims;
}

auto EquilibriumSetup::options() const -> EquilibriumOptions const&
{
    return pimpl->options;
}

auto EquilibriumSetup::Aex() const -> MatrixXdConstRef
{
    return pimpl->Aex;
}

auto EquilibriumSetup::Aen() const -> MatrixXdConstRef
{
    return pimpl->Aex.leftCols(pimpl->dims.Nn);
}

auto EquilibriumSetup::Aeq() const -> MatrixXdConstRef
{
    return pimpl->Aex.rightCols(pimpl->dims.Nq);
}

auto EquilibriumSetup::Aep() const -> MatrixXdConstRef
{
    return pimpl->Aep;
}

auto EquilibriumSetup::assembleLowerBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd
{
    return pimpl->assembleLowerBoundsVector(restrictions, state0);
}

auto EquilibriumSetup::assembleUpperBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd
{
    return pimpl->assembleUpperBoundsVector(restrictions, state0);
}

auto EquilibriumSetup::update(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> void
{
    pimpl->update(x, p, w);
}

auto EquilibriumSetup::updateGradX(VectorXlConstRef ibasicvars) -> void
{
    pimpl->updateGradX(ibasicvars);
}

auto EquilibriumSetup::updateGradP() -> void
{
    pimpl->updateGradP();
}

auto EquilibriumSetup::updateGradW() -> void
{
    pimpl->updateGradW();
}

auto EquilibriumSetup::getGibbsEnergy() -> real
{
    return pimpl->getGibbsEnergy();
}

auto EquilibriumSetup::getGibbsGradX() -> VectorXdConstRef
{
    return pimpl->getGibbsGradX();
}

auto EquilibriumSetup::getGibbsHessianX() -> MatrixXdConstRef
{
    return pimpl->getGibbsHessianX();
}

auto EquilibriumSetup::getGibbsHessianP() -> MatrixXdConstRef
{
    return pimpl->getGibbsHessianP();
}

auto EquilibriumSetup::getGibbsHessianC() -> MatrixXdConstRef
{
    return pimpl->getGibbsHessianC();
}

auto EquilibriumSetup::getConstraintResiduals() -> MatrixXdConstRef
{
    return pimpl->getConstraintResiduals();
}

auto EquilibriumSetup::getConstraintResidualsGradX() -> MatrixXdConstRef
{
    return pimpl->getConstraintResidualsGradX();
}

auto EquilibriumSetup::getConstraintResidualsGradP() -> MatrixXdConstRef
{
    return pimpl->getConstraintResidualsGradP();
}

auto EquilibriumSetup::getConstraintResidualsGradC() -> MatrixXdConstRef
{
    return pimpl->getConstraintResidualsGradC();
}

auto EquilibriumSetup::usingPartiallyExactDerivatives() -> bool
{
    return pimpl->usingPartiallyExactDerivatives();
}

auto EquilibriumSetup::usingDiagonalApproxDerivatives() -> bool
{
    return pimpl->usingDiagonalApproxDerivatives();
}

auto EquilibriumSetup::assembleChemicalPropsJacobianBegin() -> void
{
    pimpl->props.assembleFullJacobianBegin();
}

auto EquilibriumSetup::assembleChemicalPropsJacobianEnd() -> void
{
    pimpl->props.assembleFullJacobianEnd();
}

auto EquilibriumSetup::equilibriumProps() const -> EquilibriumProps const&
{
    return pimpl->props;
}

auto EquilibriumSetup::chemicalProps() const -> ChemicalProps const&
{
    return pimpl->props.chemicalProps();
}

} // namespace Reaktoro
