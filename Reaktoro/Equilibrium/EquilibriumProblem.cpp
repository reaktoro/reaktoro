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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>

namespace Reaktoro {

struct EquilibriumProblem::Impl
{
    /// The chemical system associated with this equilibrium problem
    ChemicalSystem system;

    /// The equilibrium constraints associated with this equilibrium problem
    EquilibriumConstraints constraints;

    /// The auxiliary chemical properties of the system.
    ChemicalProps props;

    Index Ne  = 0; ///< The number of chemical elements and electric charge in the chemical system.
    Index Nn  = 0; ///< The number of species in the chemical system.
    Index Npe = 0; ///< The number of equation constraints among the functional constraints.
    Index Npp = 0; ///< The number of property preservation constraints among the functional constraints.
    Index Np  = 0; ///< The number of functional constraints (Np = Npe + Npp).
    Index Nq  = 0; ///< The number of chemical potential constraints.
    Index Nir = 0; ///< The number of reactions prevented from reacting during the equilibrium calculation.
    Index Nc  = 0; ///< The number of introduced control variables (must be equal to Np)
    Index Nx  = 0; ///< The number of variables (Nx = Nn + Np + Nq).
    Index Nb  = 0; ///< The number of components (Nb = Ne + Nir).

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
      props(constraints.system())
    {
        // Initialize the auxiliary number variables
        Ne  = system.elements().size() + 1;
        Nn  = system.species().size();
        Npe = constraints.data().econstraints.size();
        Npp = constraints.data().pconstraints.size();
        Np  = Npe + Npp;
        Nq  = constraints.data().uconstraints.size();
        Nir = constraints.data().restrictions.reactions_cannot_react.size();
        Nc  = constraints.data().controls.size();
        Nx  = Nn + Np + Nq;
        Nb  = Ne + Nir;

        // Assert dimensions are proper
        error(Np != Nc, "The number of introduced control variables (", Nc, ") "
            "using method EquilibriumConstraints::control must be equal to the "
            "number of functional constraints introduced with methods "
            "EquilibriumConstraints::until (", Npe, ") and "
            "EquilibriumConstraints::preserve (", Npp, ").");

        // Initialize the auxiliary scalars, vectors and matrices
        n.resize(Nn);
        p.resize(Np);
        q.resize(Nq);
        npq.resize(Nx);
        pp0.resize(Npp);
        g.resize(Nx);
        H.resize(Nx, Nx);
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) const -> void
    {
        assert(vec.size() == Ne);
        vec[Ne - 1] = formula.charge(); // last entry in the column vector is charge of substance
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

        assert(B.rows() == Ne); // number of elements plus charge
        assert(B.cols() == Nc); // number of introduced control variables

        auto j = controls.T + controls.P; // skip columns BT and BP (if applicable), since these are zeros
        for(const auto& formula : controls.titrants)
            assembleFormulaVector(B.col(j++), formula);
    }

    /// Assemble the matrix block U in the conservation matrix C.
    auto assembleMatrixU(MatrixXdRef U) const -> void
    {
        const auto& uconstraints = constraints.data().uconstraints;

        assert(U.rows() == Ne); // number of elements plus charge
        assert(U.cols() == Nq); // number of introduced chemical potential constraints

        auto j = 0;
        for(const auto& [formula, _] : uconstraints)
            assembleFormulaVector(U.col(j++), formula);
    }

    /// Assemble the matrix block S in the conservation matrix C.
    auto assembleMatrixS(MatrixXdRef S) const -> void
    {
        assert(S.rows() == Nir);

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
        MatrixXd C = MatrixXd::Zero(Nb, Nx);

        auto A = C.topRows(Ne).leftCols(Nn);
        auto B = C.topRows(Ne).middleCols(Nn, Np);
        auto U = C.topRows(Ne).rightCols(Nq);
        auto S = C.bottomLeftCorner(Nir, Nn);

        assembleMatrixA(A);
        assembleMatrixB(B);
        assembleMatrixU(U);
        assembleMatrixS(S);

        return C;
    }

    /// Return the objective function to be minimized based on the given equilibrium constraints.
    auto objective(const ChemicalProps& props0)
    {
        // AUXILIARY REFERENCES
        const auto& controls = constraints.data().controls;
        const auto& uconstraints = constraints.data().uconstraints;
        const auto& econstraints = constraints.data().econstraints;
        const auto& pconstraints = constraints.data().pconstraints;

        // INITIALIZE TEMPERATURE AND PRESSURE IN CASE THEY ARE PRESERVED/CONSTANT
        T = props0.temperature();
        P = props0.pressure();

        // INITIALIZE THE VALUES OF PROPERTIES THAT MUST BE PRESERVED
        for(auto i = 0; i < Npp; ++i)
            pp0[i] = pconstraints[i](props0); // property value from initial chemical properties (props0)

        // DEFINE THE CHEMICAL PROPERTIES UPDATE FUNCTION
        auto update_props = [=](const VectorXr& npq)
        {
            n = npq.head(Nn);
            p = npq.segment(Nn, Np);

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
            return (n * u).sum();
        };

        // CREATE THE OBJECTIVE GRADIENT FUNCTION
        auto objective_g = [=](const VectorXr& npq)
        {
            update_props(npq); // update the chemical properties of the system with updated n, p, q

            auto gn = g.head(Nn);        // where we set the chemical potentials of the species
            auto gp = g.segment(Nn, Np); // where we set the residuals of functional constraints
            auto gq = g.tail(Nq);        // where we set the desired chemical potentials of certain substances

            auto gpe = gp.head(Npe); // where we set the residuals of the equation constraints
            auto gpp = gp.tail(Npp); // where we set the residuals of the property preservation constraints

            gn = props.chemicalPotentials(); // set the current chemical potentials of species

            for(auto i = 0; i < Nq; ++i)
                gq[i] = uconstraints[i].fn(T, P); // set the fixed chemical potentials using current T and P

            EquilibriumEquationArgs args{props, q, controls.titrants};
            for(auto i = 0; i < Npe; ++i)
                gpe[i] = econstraints[i](args); // set the residuals of the equation constraints

            for(auto i = 0; i < Npp; ++i)
                gpp[i] = pconstraints[i](props) - pp0[i]; // set the residuals of the property preservation constraints

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
        obj.f = [=](VectorXdConstRef x) {
            npq = x;
            return objective_f(npq);
            };
        obj.g = [=](VectorXdConstRef x, VectorXdRef res) { npq = x; res = objective_g(npq); };
        obj.H = [=](VectorXdConstRef x, MatrixXdRef res) { npq = x; res = objective_H(npq); };

        return obj;
    }

    auto update(const EquilibriumConstraints& constraints) -> void
    {
        error(true, "EquilibriumProblem::update not implemented yet.");
    }
};

EquilibriumProblem::EquilibriumProblem(const EquilibriumConstraints& constraints)
: pimpl(new Impl(constraints))
{}

// EquilibriumProblem::EquilibriumProblem(const EquilibriumProblem& other)
// : pimpl(new Impl(*other.pimpl))
// {}

// EquilibriumProblem::~EquilibriumProblem()
// {}

// auto EquilibriumProblem::operator=(EquilibriumProblem other) -> EquilibriumProblem&
// {
//     pimpl = std::move(other.pimpl);
//     return *this;
// }

auto EquilibriumProblem::numComponents() const -> Index
{
    return pimpl->Nb;
}

auto EquilibriumProblem::numVariables() const -> Index
{
    return pimpl->Nx;
}

auto EquilibriumProblem::conservationMatrix() const -> MatrixXd
{
    return pimpl->conservationMatrix();
}

auto EquilibriumProblem::objective(const ChemicalProps& props0) const -> EquilibriumObjective
{
    return pimpl->objective(props0);
}

auto EquilibriumProblem::update(const EquilibriumConstraints& constraints) -> void
{
    pimpl->update(constraints);
}

} // namespace Reaktoro
