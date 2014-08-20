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

#include "EquilibriumOutput.hpp"

// C++ includes
#include <iomanip>

// Eigen includes
#include <Eigen/LU>

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Equilibrium/BalanceConstraints.hpp>
#include <Reaktor/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktor/Equilibrium/EquilibriumLagrange.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>

namespace Reaktor {

auto dualChemicalPotentials(const ChemicalState& state, const EquilibriumLagrange& lagrange) -> std::tuple<Vector, Vector>
{
    // A reference to the chemical system
    const ChemicalSystem& system = lagrange.system();

    // A reference to the partitioning of the chemical species
    const Partitioning& partitioning = lagrange.partitioning();

    // The temperature of the chemical state (in units of K)
    const double T = state.temperature();

    // The pressure of the chemical state (in units of bar)
    const double P = state.pressure();

    // The product of the universal gas constant times temperature
    const double RT = universalGasConstant * T;

    // The mass-balance and charge-balance constraints
    BalanceConstraints balance(system, partitioning);

    // The mass-balance and charge-balance matrix
    const Matrix& He = balance.balanceMatrix();

    // The indices of the stable species (global indices)
    const Indices idx_stable_species = lagrange.idxStableSpecies();

    // The indices of the stable species (local indices relative to the equilibrium species)
    const Indices idx_equilibrium_stable_species = lagrange.idxEquilibriumStableSpecies();

    // The mass-balance and charge-balance matrix w.r.t. the stable species
    const Matrix Hs = cols(idx_equilibrium_stable_species, He);

    // The standard chemical potential of the species
    const Vector u0 = system.chemicalPotentials(T, P);

    // The natural log of the activities of the species
    const Vector ln_a = func(state.activities()).array().log();

    // The chemical potential of the species
    const Vector u = u0 + RT * ln_a;

    // The chemical potential of the stable species
    const Vector us = rows(idx_stable_species, u);

    // The Lagrange multipliers of the mass-balance and charge-balance constraints
    const Vector y = Hs.transpose().fullPivLu().solve(us);

    // The number of elements in the system
    const unsigned E = system.numElements();

    // The number of phases in the system
    const unsigned Np = system.numPhases();

    // The number of elements in the equilibrium partitioning
    const unsigned Ee  = partitioning.numEquilibriumElements();

    // The Lagrange multipliers of the mass-balance constraints
    const Vector y1 = y.segment(00, Ee);

    // The Lagrange multipliers of the charge-balance constraints
    const Vector y2 = y.segment(Ee, y.rows() - Ee);

    // The extended mass-balance Lagrange multipliers of the all elements (zero for non-equilibrium elements)
    Vector ye = zeros(E);
    setRows(partitioning.idxEquilibriumElements(), y1, ye);

    // The extended charge-balance Lagrange multipliers of the all phases (zero for non-equilibrium phases)
    Vector yc = zeros(Np);
    setRows(balance.idxChargeImbalancedPhases(), y2, yc);

    return std::make_tuple(ye, yc);
}

auto operator<<(std::ostream& out, const EquilibriumOutput& instance) -> std::ostream&
{
    const ChemicalState& state = instance.state;

    const EquilibriumLagrange& lagrange = instance.lagrange;

    const EquilibriumConstraints& constraints = instance.constraints;

    // A reference to the chemical system
    const ChemicalSystem& system = lagrange.system();

    // A reference to the partitioning of the chemical species
    const Partitioning& partitioning = lagrange.partitioning();

    // The mass-balance and charge-balance constraints
    BalanceConstraints balance(system, partitioning);

    const double T = state.temperature();

    const double RT = universalGasConstant * T;

    // The residuals of the mass-balance equations of the equilibrium elements
    const Vector re = func(constraints(state));

    // The number of elements in the system
    const unsigned E = system.numElements();

    // The number of phases in the system
    const unsigned Np = system.numPhases();

    // The number of elements in the equilibrium partitioning
    const unsigned Ee  = partitioning.numEquilibriumElements();

    // The Lagrange multipliers of the mass-balance and charge-balance constraints
    Vector y = lagrange.lagrangeY();

    // The Lagrange multipliers of the mass-balance constraints
    const Vector y1 = y.segment(00, Ee);

    // The Lagrange multipliers of the charge-balance constraints
    const Vector y2 = y.segment(Ee, y.rows() - Ee);

    // The extended mass-balance Lagrange multipliers of the all elements (zero for non-equilibrium elements)
    Vector ye = zeros(E);
    setRows(partitioning.idxEquilibriumElements(), y1, ye);

    // The extended charge-balance Lagrange multipliers of the all phases (zero for non-equilibrium phases)
    Vector yc = zeros(Np);
    setRows(balance.idxChargeImbalancedPhases(), y2, yc);

    // The residuals of the mass-balance equations of all elements
    Vector residuals = zeros(system.numElements());
    setRows(partitioning.idxEquilibriumElements(), re, residuals);

    // Auxiliary variables for pretty-printing
    const unsigned nfill = 5 * 25;
    const std::string bar(nfill, '=');

    const Vector n = state.composition();

    // Output the header of the element-related state
    out << bar << std::endl;
    out << std::setw(25) << std::left << "Element";
    out << std::setw(25) << std::left << "Amount";
    out << std::setw(25) << std::left << "Molality";
    out << std::setw(25) << std::left << "Residual";
    out << std::setw(25) << std::left << "Potential/Element" << std::endl;
    out << bar << std::endl;

    unsigned i = 0;
    for(std::string element : system.elements())
    {
        out << std::setw(25) << std::left  << element;
        out << std::setw(25) << std::left << state.amountElement(element)();
        out << std::setw(25) << std::left << state.molalityElement(element)();
        out << std::setw(25) << std::left << residuals[i];
        out << std::setw(25) << std::left << -ye[i]/RT;
        out << std::endl;

        ++i;
    }

    out << std::endl;

    // Output the header of the phase-related state
    out << bar << std::endl;
    out << std::setw(25) << std::left << "Phase";
    out << std::setw(25) << std::left << "Number of Species";
    out << std::setw(25) << std::left << "Amount";
    out << std::setw(25) << std::left << "Stability Index";
    out << std::setw(25) << std::left << "Potential/Charge" << std::endl;
    out << bar << std::endl;

    for(unsigned i = 0; i < system.numPhases(); ++i)
    {
        const Phase& phase = system.phase(i);

        out << std::setw(25) << std::left  << phase.name();
        out << std::setw(25) << std::left << phase.numSpecies();
        out << std::setw(25) << std::left << state.amountPhase(i)();
        out << std::setw(25) << std::left << lagrange.stabilityIndicesPhases()[i];
        out << std::setw(25) << std::left << yc[i]/RT;
        out << std::endl;
    }

    out << std::endl;

    // Output the chemical state
    out << state << std::endl;

    // Calculate pH, pe, and Eh
    const double F = 96485.3365; // the Faraday constant
    const double pH = state.acidity(state.activities());
    const double pe = yc.rows() ? - yc[0]/RT/std::log(10) : 0.0;
    const double Eh = yc.rows() ? - yc[0]/F : 0.0;

    // Calculate the ionic strength
    const double I = 0.5 * (state.composition().array() * system.speciesCharges().array().pow(2)).sum();

    // Output the pH, pe and Eh
    out << bar << std::endl;
    out << std::setw(25) << std::left << "pH";
    out << std::setw(25) << std::left << "pe";
    out << std::setw(25) << std::left << "Eh";
    out << std::setw(25) << std::left << "Ionic Strength" << std::endl;
    out << bar << std::endl;

    out << std::setw(25) << std::left << pH;
    out << std::setw(25) << std::left << pe;
    out << std::setw(25) << std::left << Eh;
    out << std::setw(25) << std::left << I << std::endl;

    return out;
}

} /* namespace Reaktor */
