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

// Optima includes
#include <Optima/Problem.hpp>

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
    VectorXr vp;   ///< The residuals of the equation constraints.
    MatrixXd Vpx;  ///< The Jacobian of the residuals of the equation constraints with respect to x = (n, q).
    MatrixXd Vpp;  ///< The Jacobian of the residuals of the equation constraints with respect to p.

    /// Construct an EquilibriumSetup::Impl object
    Impl(const EquilibriumSpecs& specs)
    : system(specs.system()), specs(specs), dims(specs), props(specs.system())
    {
        const auto Nn = dims.Nn;
        const auto Np = dims.Np;
        const auto Nq = dims.Nq;
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
        const auto Nq = dims.Nq;
        const auto Nc = dims.Nc;
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
        const auto Np = dims.Np;
        const auto Nc = dims.Nc;

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
        const auto n0 = state0.speciesAmounts().matrix();
        VectorXd be = Wn * n0;
        return be;
    }

    /// Assemble the lower bound vector `xlower` in the optimization problem where *x = (n, q)*.
    auto assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
    {
        const auto Nn = dims.Nn;
        const auto Nq = dims.Nq;
        const auto Nx = Nn + Nq;
        VectorXd xlower = constants(Nx, -inf);
        auto nlower = xlower.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.indicesSpeciesCannotDecreaseBelow()) nlower[i] = val;
        for(auto i : restrictions.indicesSpeciesCannotDecrease()) nlower[i] = n0[i]; // this comes after, in case a species cannot strictly decrease
        return xlower;
    }

    /// Assemble the upper bound vector `xupper` in the optimization problem where *x = (n, q)*.
    auto assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd
    {
        const auto Nn = dims.Nn;
        const auto Nq = dims.Nq;
        const auto Nx = Nn + Nq;
        VectorXd xupper = constants(Nx, inf);
        auto nupper = xupper.head(Nn);
        const auto n0 = state0.speciesAmounts();
        for(auto [i, val] : restrictions.indicesSpeciesCannotIncreaseAbove()) nupper[i] = val;
        for(auto i : restrictions.indicesSpeciesCannotIncrease()) nupper[i] = n0[i]; // this comes after, in case a species cannot strictly increase
        return xupper;
    }

    // /// Assemble the matrix block S in the conservation matrix C.
    // auto assembleMatrixS(MatrixXdRef S) const -> void
    // {
        // assert(S.rows() == dims.Nir);

        // const auto& inert_reactions = conditions.details().restrictions.reactions_cannot_react;

        // auto fill_matrix_row = [=](const auto& pairs, auto row)
        // {
        //     for(auto [ispecies, coeff] : pairs)
        //         row[ispecies] = coeff;
        // };

        // auto i = 0;
        // for(const auto& pairs : inert_reactions)
        //     fill_matrix_row(pairs, S.row(i++));
    // }

    // /// Return the objective function to be minimized based on the given equilibrium conditions.
    // auto objective(const ChemicalState& state0)
    // {
    //     // AUXILIARY REFERENCES
    //     const auto& uconstraints = details.uconstraints;
    //     const auto& econstraints = details.econstraints;
    //     const auto& pconstraints = details.pconstraints;

    //     const auto unknownT = specs.isTemperatureUnknown();
    //     const auto unknownP = specs.isPressureUnknown();
    //     const auto dims = this->dims;

    //     const auto tau = options.epsilon * options.logarithm_barrier_factor;

    //     // INITIALIZE TEMPERATURE AND PRESSURE IN CASE THEY ARE CONSTANT
    //     T = state0.temperature();
    //     P = state0.pressure();

    //     // // INITIALIZE THE VALUES OF PROPERTIES THAT MUST BE PRESERVED
    //     // for(auto i = 0; i < dims.Npp; ++i)
    //     //     pp0[i] = pconstraints[i].fn(state0.props()); // property value from initial chemical properties (props0)

    //     // DEFINE THE CHEMICAL PROPERTIES UPDATE FUNCTION
    //     auto update_props = [unknownT, unknownP, dims](const VectorXr& npq)
    //     {
    //         const auto n = npq.head(dims.Nn);
    //         const auto p = npq.segment(dims.Nn, dims.Np);

    //         // Update temperature and pressure if they are unknowns
    //         if(unknownT && unknownP) {
    //             T = p[0];
    //             P = p[1];
    //         } else if(unknownT) {
    //             T = p[0];
    //         } else if(unknownP) {
    //             P = p[0];
    //         }

    //         props.update(T, P, n);
    //     };

    //     // CREATE THE OBJECTIVE FUNCTION
    //     auto objective_f = [=](const VectorXr& npq)
    //     {
    //         update_props(npq);
    //         const auto& n = props.speciesAmounts();
    //         const auto& u = props.chemicalPotentials();
    //         const auto RT = universalGasConstant * T;
    //         const auto barrier = -tau * n.log().sum();
    //         return (n * u).sum()/RT + barrier;
    //     };

    //     // CREATE THE OBJECTIVE GRADIENT FUNCTION
    //     auto objective_g = [=](const VectorXr& npq)
    //     {
    //         update_props(npq); // update the chemical properties of the system with updated n, p, q

    //         const auto& n = props.speciesAmounts();
    //         const auto& u = props.chemicalPotentials();
    //         const auto RT = universalGasConstant * T;

    //         auto gn = g.head(dims.Nn);        // where we set the chemical potentials of the species
    //         auto gp = g.segment(dims.Nn, dims.Np); // where we set the residuals of functional conditions
    //         auto gq = g.tail(dims.Nq);        // where we set the desired chemical potentials of certain substances

    //         auto gpe = gp.head(dims.Npe); // where we set the residuals of the equation conditions
    //         auto gpp = gp.tail(dims.Npp); // where we set the residuals of the property preservation conditions

    //         gn = u/RT - tau/n; // set the current chemical potentials of species

    //         for(auto i = 0; i < dims.Nq; ++i)
    //             gq[i] = uconstraints[i].fn(T, P); // set the fixed chemical potentials using current T and P

    //         for(auto i = 0; i < dims.Npe; ++i)
    //             gpe[i] = econstraints[i].fn(props); // set the residuals of the equation conditions

    //         for(auto i = 0; i < dims.Npp; ++i)
    //             gpp[i] = pconstraints[i].fn(props) - pp0[i]; // set the residuals of the property preservation conditions

    //         return g;
    //     };

    //     // CREATE THE OBJECTIVE HESSIAN FUNCTION
    //     auto objective_H = [=](VectorXr& npq)
    //     {
    //         H = autodiff::jacobian(objective_g, autodiff::wrt(npq), autodiff::at(npq));
    //         return H;
    //     };

    //     // CONSTRUCT THE OBJECTIVE FUNCTION, ITS GRADIENT AND HESSIAN
    //     EquilibriumObjective obj;
    //     obj.f = [=](VectorXdConstRef x) { npq = x; return objective_f(npq); };
    //     obj.g = [=](VectorXdConstRef x, VectorXdRef res) { npq = x; res = objective_g(npq); };
    //     obj.H = [=](VectorXdConstRef x, MatrixXdRef res) { npq = x; res = objective_H(npq); };

    //     return obj;
    // }

    /// Update the chemical properties of the system.
    auto updateChemicalProps(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> void
    {
        //======================================================================
        // Update temperature (T) and pressure (P)
        //----------------------------------------------------------------------
        // If temperature or pressure are unknowns, get their current values
        // from vector p. If not, they have been given in the params object.
        //======================================================================
        const auto unknownT = specs.isTemperatureUnknown();
        const auto unknownP = specs.isPressureUnknown();

        // If  either from valuesif they are unknowns
        if(unknownT && unknownP)
        {
            assert(p.size() >= 2);
            T = p[0];
            P = p[1];
        }
        else if(unknownT)
        {
            assert(p.size() >= 1);
            T = p[0];
            P = params.get("P").value();
        }
        else if(unknownP)
        {
            assert(p.size() >= 1);
            T = params.get("T").value();
            P = p[0];
        }
        else
        {
            T = params.get("T").value();
            P = params.get("P").value();
        }

        assert(T > 0.0); // check if proper value for T was given in p or params
        assert(P > 0.0); // check if proper value for P was given in p or params

        //======================================================================
        // Update chemical properties of the chemical system
        //======================================================================
        const auto n = x.head(dims.Nn);

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

        const auto Nn = dims.Nn;
        const auto Nq = dims.Nq;

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
            gq[i] = uconstraints[i].fn(props, params);

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

    auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef
    {
        updateChemicalProps(x, p, params);

        const auto Np = dims.Np;

        const auto& econstraints = specs.constraintsEquationType();

        for(auto i = 0; i < Np; ++i)
            vp[i] = econstraints[i].fn(props, params);

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

} // namespace Reaktoro
