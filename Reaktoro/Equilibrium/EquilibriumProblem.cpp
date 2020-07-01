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

#include "EquilibriumProblem.hpp"

// Optima includes
#include <Optima/Problem.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>

namespace Reaktoro {

struct EquilibriumProblem::Impl
{
    /// The chemical system associated with this equilibrium problem
    ChemicalSystem system;

    /// The equilibrium constraints associated with this equilibrium problem
    EquilibriumConstraints constraints;

    /// The auxiliary chemical properties of the system.
    ChemicalProps props;

    /// The dimensions of the variables in the equilibrium problem.
    EquilibriumDims dims;

    real T = 0.0;  ///< The current temperature of the chemical system (in K).
    real P = 0.0;  ///< The current pressure of the chemical system (in Pa).
    ArrayXr n;     ///< The current amounts of the species (in mol).
    ArrayXr p;     ///< The current values of the introduced control variables for the functional constraints.
    ArrayXr q;     ///< The current amounts of the introduced substances whose chemical potentials are constrained.
    VectorXr npq;  ///< The auxiliary vector (n, p, q).
    ArrayXr pp0;   ///< The initial values of the preserved properties.
    VectorXr g;    ///< The gradient vector of the objective function.
    MatrixXd H;    ///< The Hessian matrix of the objective function.

    /// Construct an EquilibriumProblem::Impl object
    Impl(const EquilibriumConstraints& constraints)
    : system(constraints.system()),
      constraints(constraints),
      props(constraints.system()),
      dims(constraints)
    {
        // Initialize the auxiliary scalars, vectors and matrices
        n.resize(dims.Nn);
        p.resize(dims.Np);
        q.resize(dims.Nq);
        npq.resize(dims.Nx);
        pp0.resize(dims.Npp);
        g.resize(dims.Nx);
        H.resize(dims.Nx, dims.Nx);
    }

    /// Update the equilibrium constraints for the next chemical equilibrium calculation.
    auto update(const EquilibriumConstraints& _constraints) -> void
    {
        error(!_constraints.data().locked, "EquilibriumProblem::update "
            "accepts only a locked EquilibriumConstraints object.");
        constraints = _constraints;
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) const -> void
    {
        assert(vec.size() == dims.Ne);
        vec[dims.Ne - 1] = formula.charge(); // last entry in the column vector is charge of substance
        for(const auto& [element, coeff] : formula.elements()) {
            const auto ielem = system.elements().index(element);
            vec[ielem] = coeff;
        }
    }

    /// Assemble the matrix block A in the conservation matrix C.
    auto assembleMatrixA(MatrixXdRef A) const -> void
    {
        A = system.formulaMatrix();
    }

    /// Assemble the matrix block B = [BT BP Bq] in the conservation matrix C.
    auto assembleMatrixB(MatrixXdRef B) const -> void
    {
        const auto& controls = constraints.data().controls;

        assert(B.rows() == dims.Ne); // number of elements plus charge
        assert(B.cols() == dims.Np); // number of explicitly introduced control variables

        auto j = controls.T + controls.P; // skip columns BT and BP (if applicable), since these are zeros
        for(const auto& formula : controls.titrants)
            assembleFormulaVector(B.col(j++), formula);
    }

    /// Assemble the matrix block U in the conservation matrix C.
    auto assembleMatrixU(MatrixXdRef U) const -> void
    {
        const auto& uconstraints = constraints.data().uconstraints;

        assert(U.rows() == dims.Ne); // number of elements plus charge
        assert(U.cols() == dims.Nq); // number of introduced chemical potential constraints

        auto j = 0;
        for(const auto& [formula, _] : uconstraints)
            assembleFormulaVector(U.col(j++), formula);
    }

    /// Assemble the matrix block S in the conservation matrix C.
    auto assembleMatrixS(MatrixXdRef S) const -> void
    {
        assert(S.rows() == dims.Nir);

        const auto& inert_reactions = constraints.data().restrictions.reactions_cannot_react;

        auto fill_matrix_row = [=](const auto& pairs, auto row)
        {
            for(auto [ispecies, coeff] : pairs)
                row[ispecies] = coeff;
        };

        auto i = 0;
        for(const auto& pairs : inert_reactions)
            fill_matrix_row(pairs, S.row(i++));
    }

    /// Assemble the conservation matrix based on the given equilibrium constraints.
    /// The conservation matrix is:
    ///
    /// C = [ A B U ]
    ///     [ S 0 0 ]
    ///
    /// where A is the formula matrix of the species with respect to elements
    /// and charge; B = [BT BP Bq], with BT and BP being zero column vectors
    /// and Bq the formula matrix of the introduced titrants whose amounts are
    /// controlled to attain imposed equilibrium constraints; U is the formula
    /// matrix of the substances with fixed chemical potentials; and S is the
    /// stoichiometric matrix of reactions that cannot progress during the
    /// equilibrium calculation (inert reactions).
    auto conservationMatrix() const -> MatrixXd
    {
        MatrixXd C = MatrixXd::Zero(dims.Nc, dims.Nx);

        auto A = C.topRows(dims.Ne).leftCols(dims.Nn);
        auto B = C.topRows(dims.Ne).middleCols(dims.Nn, dims.Np);
        auto U = C.topRows(dims.Ne).rightCols(dims.Nq);
        auto S = C.bottomLeftCorner(dims.Nir, dims.Nn);

        assembleMatrixA(A);
        assembleMatrixB(B);
        assembleMatrixU(U);
        assembleMatrixS(S);

        return C;
    }

    /// Return the objective function to be minimized based on the given equilibrium constraints.
    auto objective(const ChemicalState& state0)
    {
        // AUXILIARY REFERENCES
        const auto& controls = constraints.data().controls;
        const auto& uconstraints = constraints.data().uconstraints;
        const auto& econstraints = constraints.data().econstraints;
        const auto& pconstraints = constraints.data().pconstraints;

        // INITIALIZE TEMPERATURE AND PRESSURE IN CASE THEY ARE PRESERVED/CONSTANT
        T = state0.temperature();
        P = state0.pressure();

        // INITIALIZE THE VALUES OF PROPERTIES THAT MUST BE PRESERVED
        for(auto i = 0; i < dims.Npp; ++i)
            pp0[i] = pconstraints[i].fn(state0.props()); // property value from initial chemical properties (props0)

        // DEFINE THE CHEMICAL PROPERTIES UPDATE FUNCTION
        auto update_props = [=](const VectorXr& npq)
        {
            n = npq.head(dims.Nn);
            p = npq.segment(dims.Nn, dims.Np);

            // Update temperature and pressure if they are variable
            if(controls.T && controls.P) {
                T = p[0];
                P = p[1];
            } else if(controls.T) {
                T = p[0];
            } else if(controls.P) {
                P = p[0];
            }

            props.update(T, P, n);
        };

        // CREATE THE OBJECTIVE FUNCTION
        auto objective_f = [=](const VectorXr& npq)
        {
            update_props(npq);
            const auto& n = props.speciesAmounts();
            const auto& u = props.chemicalPotentials();
            const auto RT = universalGasConstant * T;
            return (n * u).sum()/RT;
        };

        // CREATE THE OBJECTIVE GRADIENT FUNCTION
        auto objective_g = [=](const VectorXr& npq)
        {
            update_props(npq); // update the chemical properties of the system with updated n, p, q

            auto gn = g.head(dims.Nn);        // where we set the chemical potentials of the species
            auto gp = g.segment(dims.Nn, dims.Np); // where we set the residuals of functional constraints
            auto gq = g.tail(dims.Nq);        // where we set the desired chemical potentials of certain substances

            auto gpe = gp.head(dims.Npe); // where we set the residuals of the equation constraints
            auto gpp = gp.tail(dims.Npp); // where we set the residuals of the property preservation constraints

            const auto RT = universalGasConstant * T;

            gn = props.chemicalPotentials(); // set the current chemical potentials of species
            gn /= RT;

            for(auto i = 0; i < dims.Nq; ++i)
                gq[i] = uconstraints[i].fn(T, P); // set the fixed chemical potentials using current T and P

            EquilibriumEquationArgs args{props, q, controls.titrants};
            for(auto i = 0; i < dims.Npe; ++i)
                gpe[i] = econstraints[i].fn(args); // set the residuals of the equation constraints

            for(auto i = 0; i < dims.Npp; ++i)
                gpp[i] = pconstraints[i].fn(props) - pp0[i]; // set the residuals of the property preservation constraints

            return g;
        };

        // CREATE THE OBJECTIVE HESSIAN FUNCTION
        auto objective_H = [=](VectorXr& npq)
        {
            H = autodiff::jacobian(objective_g, autodiff::wrt(npq), autodiff::at(npq));
            return H;
        };

        // CONSTRUCT THE OBJECTIVE FUNCTION, ITS GRADIENT AND HESSIAN
        EquilibriumObjective obj;
        obj.f = [=](VectorXdConstRef x) { npq = x; return objective_f(npq); };
        obj.g = [=](VectorXdConstRef x, VectorXdRef res) { npq = x; res = objective_g(npq); };
        obj.H = [=](VectorXdConstRef x, MatrixXdRef res) { npq = x; res = objective_H(npq); };

        return obj;
    }

    /// Set the lower bounds of the species amounts based on the given equilibrium constraints.
    auto xlower(const ChemicalState& state0, ArrayXdRef res) const -> void
    {
        auto nlower = res.head(dims.Nn);
        const auto& n = state0.speciesAmounts();
        const auto& restrictions = constraints.data().restrictions;
        for(auto i : restrictions.species_cannot_decrease) nlower[i] = n[i];
        for(auto i : restrictions.species_cannot_react) nlower[i] = n[i];
    }

    /// Set the upper bounds of the species amounts based on the given equilibrium constraints.
    auto xupper(const ChemicalState& state0, ArrayXdRef res) const -> void
    {
        auto nupper = res.head(dims.Nn);
        const auto& n = state0.speciesAmounts();
        const auto& restrictions = constraints.data().restrictions;
        for(auto i : restrictions.species_cannot_increase) nupper[i] = n[i];
        for(auto i : restrictions.species_cannot_react) nupper[i] = n[i];
    }
};

EquilibriumProblem::EquilibriumProblem(const EquilibriumConstraints& constraints)
: pimpl(new Impl(constraints))
{}

EquilibriumProblem::EquilibriumProblem(const EquilibriumProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumProblem::~EquilibriumProblem()
{}

auto EquilibriumProblem::operator=(EquilibriumProblem other) -> EquilibriumProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumProblem::update(const EquilibriumConstraints& constraints) -> void
{
    pimpl->update(constraints);
}

auto EquilibriumProblem::dims() const -> const EquilibriumDims&
{
    return pimpl->dims;
}

auto EquilibriumProblem::conservationMatrix() const -> MatrixXd
{
    return pimpl->conservationMatrix();
}

auto EquilibriumProblem::objective(const ChemicalState& state0) const -> EquilibriumObjective
{
    return pimpl->objective(state0);
}

auto EquilibriumProblem::xlower(const ChemicalState& state0, ArrayXdRef res) const -> void
{
    pimpl->xlower(state0, res);
}

auto EquilibriumProblem::xupper(const ChemicalState& state0, ArrayXdRef res) const -> void
{
    pimpl->xupper(state0, res);
}

} // namespace Reaktoro
