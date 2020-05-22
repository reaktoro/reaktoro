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

    /// Construct an EquilibriumProblem::Impl object
    Impl(const EquilibriumConstraints& constraints)
    : system(constraints.system()), constraints(constraints)
    {
    }

    /// Return the number of chemical potential constraints.
    auto numChemicalPotentialConstraints() const -> Index
    {
        return constraints.data().uconstraints.size();
    }

    /// Return the number of components associated with given equilibrium constraints.
    auto numComponents() const -> Index
    {
        const auto& restrictions = constraints.data().restrictions;

        const auto num_elements = system.elements().size();
        const auto num_charge = 1;
        const auto num_inert_reactions = restrictions.reactions_cannot_react.size();

        return num_elements + num_charge + num_inert_reactions;
    }

    /// Return the total number of variables associated with given equilibrium constraints.
    auto numVariables() const -> Index
    {
        const auto& uconstraints = constraints.data().uconstraints;
        const auto& controls = constraints.data().controls;

        const auto num_species = system.species().size();
        const auto num_uconstraints = uconstraints.size();
        const auto num_controls = controls.size();

        return num_species + num_uconstraints + num_controls;
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) const -> void
    {
        const auto num_elements = system.elements().size();
        assert(vec.size() == num_elements + 1);
        vec[num_elements] = formula.charge(); // last entry in the column vector is charge of substance
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

    /// Assemble the matrix block U in the conservation matrix C.
    auto assembleMatrixU(MatrixXdRef U) const -> void
    {
        const auto& uconstraints = constraints.data().uconstraints;

        const auto num_elements = system.elements().size();
        const auto num_substances = uconstraints.size();

        assert(U.rows() == 1 + num_elements); // number of elements plus charge
        assert(U.cols() == num_substances);   // number of introduced chemical potential constraints

        auto j = 0;
        for(const auto& [formula, _] : uconstraints)
            assembleFormulaVector(U.col(j++), formula);
    }

    /// Assemble the matrix block B = [BT BP Bq] in the conservation matrix C.
    auto assembleMatrixB(MatrixXdRef B) const -> void
    {
        const auto& controls = constraints.data().controls;

        const auto num_elements = system.elements().size();
        const auto num_controls = controls.size();

        assert(B.rows() == 1 + num_elements); // number of elements plus charge
        assert(B.cols() == num_controls);     // number of introduced control variables

        auto j = controls.T + controls.P; // skip columns BT and BP (if applicable), since these are zeros
        for(const auto& formula : controls.titrants)
            assembleFormulaVector(B.col(j++), formula);
    }

    /// Assemble the matrix block S in the conservation matrix C.
    auto assembleMatrixS(MatrixXdRef S) const -> void
    {
        const auto& inert_reactions = constraints.data().restrictions.reactions_cannot_react;

        assert(S.rows() == inert_reactions.size());

        auto fill_matrix_row = [&](const auto& pairs, auto row)
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
    /// C = [ A U B ]
    ///     [ S 0 0 ]
    ///
    /// where A is the formula matrix of the species with respect to elements
    /// and charge; U is the formula matrix of the substances with fixed
    /// chemical potentials; B = [BT BP Bq], with BT and BP being zero column
    /// vectors and Bq the formula matrix of the introduced titrants whose
    /// amounts are controlled to attain imposed equilibrium constraints; and S
    /// is the stoichiometric matrix of reactions that cannot progress during
    /// the equilibrium calculation (inert reactions).
    auto conservationMatrix() const -> MatrixXd
    {
        const auto& restrictions = constraints.data().restrictions;
        const auto& uconstraints = constraints.data().uconstraints;
        const auto& controls = constraints.data().controls;

        const auto num_elements = system.elements().size();
        const auto num_species = system.species().size();
        const auto num_charge = 1;
        const auto num_inert_reactions = restrictions.reactions_cannot_react.size();
        const auto num_fixed_chemical_potentials = uconstraints.size();
        const auto num_controls = controls.size();

        const auto num_rows = num_elements + num_charge + num_inert_reactions;
        const auto num_cols = num_species + num_fixed_chemical_potentials + num_controls;

        MatrixXd C = MatrixXd::Zero(num_rows, num_cols);

        auto A = C.topRows(num_elements + num_charge).leftCols(num_species);
        auto U = C.topRows(num_elements + num_charge).middleCols(num_species, num_fixed_chemical_potentials);
        auto B = C.topRows(num_elements + num_charge).rightCols(num_controls);
        auto S = C.bottomLeftCorner(num_inert_reactions, num_species);

        assembleMatrixA(A);
        assembleMatrixU(U);
        assembleMatrixB(B);
        assembleMatrixS(S);

        return C;
    }

    /// Assemble the objective function to be minimized based on the given equilibrium constraints.
    auto objective() const -> EquilibriumObjective
    {
        EquilibriumObjective obj;

        const auto& uconstraints = constraints.data().uconstraints;
        const auto& econstraints = constraints.data().econstraints;
        const auto& pconstraints = constraints.data().pconstraints;

        const auto Nx = system.species().size();
        const auto Nq = constraints.numChemicalPotentialConstraints();
        const auto Np = constraints.numControlVariables();

        const auto Nec = constraints.numEquationConstraints();
        const auto Npc = constraints.numPropertyPreservationConstraints();

        error(Np != Nec + Npc,
            "Could not assemble the objective function for the "
            "chemical equilibrium calculation with given constraints. "
            "The number of control variables, ", Np, " must match the sum of "
            "the number of equation constraints, ", Nec, " and "
            "the number of property preservation constraints, ", Npc, ".");

        obj.f = [](const ChemicalProps& props)
        {
            const auto& n = props.speciesAmounts();
            const auto& u = props.chemicalPotentials();
            return (n * u).sum();
        };

        obj.g = [=](const ChemicalProps& props, VectorXrRef res)
        {
            // assert(res.size() == Nx + Nq + Np);
            assert(res.size() == Nx);

            const auto& u = props.chemicalPotentials();

            auto gx = res.head(Nx);
            // auto gq = res.segment(Nx, Nq);
            // auto gp = res.tail(Np);

            // auto gec = gp.head(Nec);
            // auto gpc = gp.tail(Npc);

            // const auto T = props.temperature();
            // const auto P = props.pressure();

            gx = u;

            // for(auto i = 0; i < Nq; ++i)
            //     gq[i] = uconstraints[i].fn(T, P);

            // for(auto i = 0; i < Nec; ++i)
            //     gec[i] = econstraints[i]({props, });

            // gp = u;

            // assert(res.size() == numSpecies() + num)
            // const auto& n = props.speciesAmounts();
            // return (n * u).sum();
        };

        obj.H = [=](const ChemicalProps& props, MatrixXdRef res)
        {
            // assert(res.size() == Nx + Nq + Np);

            // const auto& u = props.chemicalPotentials();

            // auto gx = res.head(Nx);
            // auto gq = res.segment(Nx, Nq);
            // auto gp = res.tail(Np);

            // auto gec = gp.head(Nec);
            // auto gpc = gp.tail(Npc);

            // const auto T = props.temperature();
            // const auto P = props.pressure();

            // for(auto )
            // gx = u;

            // for(auto i = 0; i < Nq; ++i)
            //     gq[i] = uconstraints[i].fn(T, P);

            // for(auto i = 0; i < Nec; ++i)
            //     gec[i] = econstraints[i]({props, });

            // gp = u;

            // assert(res.size() == numSpecies() + num)
            // const auto& n = props.speciesAmounts();
            // return (n * u).sum();
        };
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

auto EquilibriumProblem::conservationMatrix() const -> MatrixXd
{
    return pimpl->conservationMatrix();
}

auto EquilibriumProblem::objective() const -> EquilibriumObjective
{

}

} // namespace Reaktoro
