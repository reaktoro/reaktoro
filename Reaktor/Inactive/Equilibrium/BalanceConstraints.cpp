/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "BalanceConstraints.hpp"

// Eigen includes
#include <Eigen/Dense>

// Reaktor includes
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>

namespace Reaktor {

class BalanceConstraints::Impl
{
private:
    /// The indices of the equilibrium species
    Indices idx_equilibrium_species$;

    /// The indices of the equilibrium species
    Indices idx_kinetic_species$;

    /// The indices of the equilibrium species
    Indices idx_inert_species$;

    /// The indices of the phases containing equilibrium species whose charge balance condition is not implicitly imposed
    Indices idx_charge_imbalanced_phases;

    /// The formula matrix of the equilibrium species
    Matrix We$;

    /// The combined mass and charge balance matrix of the equilibrium species
    Matrix He$;

    /// The charge balance matrix of the equilibrium, kinetic and inert species
    Matrix Ce$, Ck$, Ci$;

    /// The molar abundance of the equilibrium, kinetic and inert species
    Vector ne$, nk$, ni$;

    /// The molar abundance of the equilibrium elements
    Vector be$;

    /// The result of the equilibrium constraint function evaluation
    VectorResult he$;

    /// The number of mass-balance constraints that need to be imposed in the equilibrium system
    unsigned num_mass_constraints$;

    /// The number of charge-balance constraints that need to be imposed in the equilibrium system
    unsigned num_charge_constraints$;

public:
    Impl(const ChemicalSystem& system)
    : Impl(system, Partitioning(system))
    {}

    Impl(const ChemicalSystem& system, const Partitioning& partitioning)
    {
        // Define auxiliary variables
        const unsigned N  = system.numSpecies();
        const unsigned Np = system.numPhases();
        const unsigned Ne = partitioning.numEquilibriumSpecies();
        const unsigned Ee = partitioning.numEquilibriumElements();

        // Initialise the indices of the equilibrium, kinetic and inert species
        idx_equilibrium_species$ = partitioning.idxEquilibriumSpecies();
        idx_kinetic_species$ = partitioning.idxKineticSpecies();
        idx_inert_species$ = partitioning.idxInertSpecies();

        // Initialise the formula matrix of the equilibrium species
        We$ = partitioning.equilibriumFormulaMatrix(system);

        // Assemble the matrix of electrical charges for every phase
        Matrix Z = zeros(Np, N);
        for(unsigned i = 0; i < Np; ++i)
            for(unsigned j : system.idxSpeciesInPhase(i))
                Z(i, j) = system.species(j).charge();

        // Determine the matrices of electrical charges of the phases associated with the equilibrium, kinetic and inert species
        const Matrix Ze = partitioning.equilibriumCols(Z);
        const Matrix Zk = partitioning.kineticCols(Z);
        const Matrix Zi = partitioning.inertCols(Z);

        // Determine which phases are not implicitly charge balanced by checking the linearly dependent rows of Ze w.r.t. We
        Matrix Xe = zeros(Ee + 1, Ne);
        Xe.topRows(Ee) = We$;

        for(unsigned i = 0; i < Np; ++i)
        {
            // Set the last row of Xe to the i-th row of matrix Ze
            Xe.bottomRows(1) = Ze.row(i);

            // Check if the last row is now
            if(Xe.fullPivLu().rank() == Xe.rows())
                idx_charge_imbalanced_phases.push_back(i);
        }

        // Initialise the charge balance matrix of the equilibrium, kinetic and inert species
        Ce$ = rows(idx_charge_imbalanced_phases, Ze);
        Ck$ = rows(idx_charge_imbalanced_phases, Zk);
        Ci$ = rows(idx_charge_imbalanced_phases, Zi);

        // Initialise the combined mass and charge balance matrix
        He$ = zeros(We$.rows() + Ce$.rows(), Ne);
        He$.topRows(We$.rows()) = We$;
        He$.bottomRows(Ce$.rows()) = Ce$;

        // Initialise the result of the equilibrium constraint function evaluation
        func(he$).resize(We$.rows() + Ce$.rows());
        grad(he$) = He$;

        // Initialise the number of mass and charge-balance constraints
        num_mass_constraints$ = We$.rows();
        num_charge_constraints$ = Ce$.rows();
    }

    auto setMassBalance(const Vector& be) -> void
    {
        be$ = be;
    }

    auto numConstraints() const -> unsigned
    {
        return num_mass_constraints$ + num_charge_constraints$;
    }

    auto numMassBalanceConstraints() const -> unsigned
    {
        return num_mass_constraints$;
    }

    auto numChargeBalanceConstraints() const -> unsigned
    {
        return num_charge_constraints$;
    }

    auto idxChargeImbalancedPhases() const -> const Indices&
    {
        return idx_charge_imbalanced_phases;
    }

    auto balanceMatrix() const  -> const Matrix&
    {
        return He$;
    }

    auto operator()(const ChemicalState& state) -> VectorResult
    {
        ne$ = rows(idx_equilibrium_species$, state.composition());
        nk$ = rows(idx_kinetic_species$, state.composition());
        ni$ = rows(idx_inert_species$, state.composition());

        func(he$).topRows(We$.rows()) = We$*ne$ - be$;
        func(he$).bottomRows(Ce$.rows()) = Ce$*ne$ - Ck$*nk$ - Ci$*ni$;

        return he$;
    }
};

BalanceConstraints::BalanceConstraints(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

BalanceConstraints::BalanceConstraints(const ChemicalSystem& system, const Partitioning& partitioning)
: pimpl(new Impl(system, partitioning))
{}

BalanceConstraints::BalanceConstraints(const BalanceConstraints& other)
: pimpl(new Impl(*other.pimpl))
{}

BalanceConstraints::~BalanceConstraints()
{}

auto BalanceConstraints::operator=(BalanceConstraints other) -> BalanceConstraints&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto BalanceConstraints::setMassBalance(const Vector& be) -> void
{
    return pimpl->setMassBalance(be);
}

auto BalanceConstraints::numConstraints() const -> unsigned
{
    return pimpl->numConstraints();
}

auto BalanceConstraints::numMassBalanceConstraints() const -> unsigned
{
    return pimpl->numMassBalanceConstraints();
}

auto BalanceConstraints::numChargeBalanceConstraints() const -> unsigned
{
    return pimpl->numChargeBalanceConstraints();
}

auto BalanceConstraints::idxChargeImbalancedPhases() const -> const Indices&
{
    return pimpl->idxChargeImbalancedPhases();
}

auto BalanceConstraints::balanceMatrix() const  -> const Matrix&
{
    return pimpl->balanceMatrix();
}

auto BalanceConstraints::operator()(const ChemicalState& state) -> VectorResult
{
    return pimpl->operator()(state);
}

auto BalanceConstraints::constraints(const Vector& be) const -> EquilibriumConstraints
{
    BalanceConstraints balance(*this);
    balance.setMassBalance(be);
    return EquilibriumConstraints(balance, balance.numConstraints());
}

BalanceConstraints::operator EquilibriumConstraints() const
{
    EquilibriumConstraints constraints(*this, numConstraints());
    return constraints;
}

} // namespace Reaktor
