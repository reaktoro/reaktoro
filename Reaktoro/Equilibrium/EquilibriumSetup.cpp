// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include "EquilibriumSetup.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

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

    VectorXr x;    ///< The auxiliary vector x = (n, q).
    VectorXr n;    ///< The current amounts of the species (in mol).
    VectorXr q;    ///< The current values of the introduced *q* control variables (amounts of implicit titrants whose chemical potentials are constrained).
    VectorXr p;    ///< The current values of the introduced *p* control variables (temperature, pressure, amounts of explicit titrants).
    VectorXr w;    ///< The current values of the input variables for the equilibrium calculation (e.g, T, P, pH, pE, V, etc.).

    VectorXr F;    ///< The auxiliary vector F = (gx, vp) containing first-order optimality conditions (gx) and residuals of equation constraints (vp)

    real f;        ///< The objective function value (Gibbs energy of the system).

    VectorXd gx;   ///< The gradient of the objective function with respect to x = (n, q).
    VectorXd vp;   ///< The residual vector of the equilibrium equation constraints (of size p).

    MatrixXd Hxx;  ///< The Jacobian of the objective gradient function with respect to x = (n, q).
    MatrixXd Hxp;  ///< The Jacobian of the objective gradient function with respect to p.
    MatrixXd Hxc;  ///< The Jacobian of the objective gradient function with respect to c = (w, b).

    MatrixXd Vpx;  ///< The Jacobian of vp with respect to x = (n, q).
    MatrixXd Vpp;  ///< The Jacobian of vp with respect to p.
    MatrixXd Vpc;  ///< The Jacobian of vp with respect to c = (w, b).

    VectorXl isbasicvar; /// The bitmap that indicates which variables in x = (n, q) are currently basic variables.

    Indices ipps;  ///< The indices of the pure phase species (i.e., species composing single-phase species, whose chemical potentials do not depend on composition)

    // -------------------------------------------- //
    // ------ CONVENIENT AUXILIARY VARIABLES ------ //
    // -------------------------------------------- //

    const Index Nb; ///< The number of conservative components
    const Index Nn; ///< The number of chemical species
    const Index Np; ///< The number of p control variables
    const Index Nq; ///< The number of q control variables
    const Index Nx; ///< The number of x variables with x = (n, q)
    const Index Nw; ///< The number of input variables w in the chemical equilibrium problem
    const Index Nc; ///< The number of input variables in c = (w, b)

    /// Construct an EquilibriumSetup::Impl object
    Impl(const EquilibriumSpecs& specs)
    : system(specs.system()), specs(specs), dims(specs), props(specs),
      Nb(dims.Nb), Nn(dims.Nn), Np(dims.Np), Nq(dims.Nq), Nx(dims.Nx), Nw(dims.Nw), Nc(dims.Nw + dims.Nb)
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

        isbasicvar.resize(Nx);

        auto offset = 0;
        for(const auto& phase : system.phases())
        {
            const auto size = phase.species().size();
            if(size == 1)
                ipps.push_back(offset);
            offset += size;
        }
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) const -> void
    {
        assert(vec.size() == Nb);
        vec.fill(0.0);
        vec[Nb - 1] = formula.charge(); // last entry in the column vector is charge of substance
        for(const auto& [element, coeff] : formula.elements()) {
            const auto ielem = system.elements().index(element);
            vec[ielem] = coeff;
        }
    }

    /// Assemble the coefficient matrix `Aex` in optimization problem.
    auto assembleMatrixAex() const -> MatrixXd
    {
        MatrixXd Aex = zeros(Nb, Nx);

        auto Wn = Aex.topLeftCorner(Nb, Nn);  // the formula matrix of the species
        auto Wq = Aex.topRightCorner(Nb, Nq); // the formula matrix of the implicit titrants

        Wn = system.formulaMatrix();

        auto j = 0;
        for(const auto& formula : specs.titrantsImplicit())
            assembleFormulaVector(Wq.col(j++), formula);

        return Aex;
    }

    /// Assemble the coefficient matrix `Aep` in optimization problem.
    auto assembleMatrixAep() const -> MatrixXd
    {
        MatrixXd Aep = zeros(Nb, Np);

        auto Wp = Aep.topRows(Nb); // the formula matrix of temperature, pressure and explicit titrants

        auto j = specs.isTemperatureUnknown() + specs.isPressureUnknown(); // skip columns corresponding to T and P in p (if applicable), since these are zeros
        for(const auto& formula : specs.titrantsExplicit())
            assembleFormulaVector(Wp.col(j++), formula);

        return Aep;
    }

    auto assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
    {
        VectorXd xlower = constants(Nx, -inf);
        auto nlower = xlower.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.speciesCannotDecreaseBelow()) nlower[i] = val;
        for(auto i : restrictions.speciesCannotDecrease()) nlower[i] = n0[i]; // this comes after, in case a species cannot strictly decrease
        for(auto& val : nlower) val = std::max(val, options.epsilon); // ensure the upper bounds of the species amounts are not below the minimum amount value given in EquilibriumOptions::epsilon. TODO: Issue a warning when lower/upper bound of a species amount is changed to EquilibriumOptions::epsilon.
        return xlower;
    }

    auto assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
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

        props.update(n, p, w);

        updateGibbsEnergy();
        updateF();
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
        const auto& uconstraints = specs.constraintsChemicalPotentialType();
        const auto& econstraints = specs.constraintsEquationType();

        const auto& cprops = props.chemicalProps();

        const auto& T  = cprops.temperature();
        const auto& n  = cprops.speciesAmounts();
        const auto& mu = cprops.chemicalPotentials();

        const auto RT  = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;

        auto sn = F.head(Nn);        // the segment in F where we set the chemical potentials of the species
        auto sq = F.segment(Nn, Nq); // the segment in F where we set the desired chemical potentials of some substances
        auto sp = F.tail(Np);        // the segment in F where we set the residuals of the equation constraints

        sn = mu/RT; // set the current chemical potentials of species (normalized by RT)
        sn(ipps).array() -= tau/n(ipps); // add log barrier contribution to pure phase species

        for(auto i = 0; i < Nq; ++i)
            sq[i] = uconstraints[i].fn(cprops, w)/RT;

        for(auto i = 0; i < Np; ++i)
            sp[i] = econstraints[i].fn(cprops, w);
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
        const auto& cprops = props.chemicalProps();
        const auto& T = cprops.temperature();
        const auto& n = cprops.speciesAmounts();
        const auto& u = cprops.chemicalPotentials();
        const auto RT = universalGasConstant * T;
        const auto tau = options.epsilon * options.logarithm_barrier_factor;
        const auto barrier = -tau * n(ipps).log().sum();
        f = (n * u).sum()/RT + barrier; // the current Gibbs energy of the system (normalized by RT)
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

auto EquilibriumSetup::assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
{
    return pimpl->assembleLowerBoundsVector(restrictions, state0);
}

auto EquilibriumSetup::assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
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

auto EquilibriumSetup::equilibriumProps() const -> const EquilibriumProps&
{
    return pimpl->props;
}

auto EquilibriumSetup::chemicalProps() const -> const ChemicalProps&
{
    return pimpl->props.chemicalProps();
}

} // namespace Reaktoro
