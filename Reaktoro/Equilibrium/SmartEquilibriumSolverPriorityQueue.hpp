// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#pragma once

// C++ includes
#include <deque>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolverBase.hpp>

namespace Reaktoro {
    
/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolverPriorityQueue : public SmartEquilibriumSolverBase
{
public:

    /// Construct an SmartEquilibriumSolverPriorityQueue instance with given chemical system.
    SmartEquilibriumSolverPriorityQueue(const ChemicalSystem& system);

    /// Construct an SmartEquilibriumSolverPriorityQueue instance with given partition of the chemical system.
    SmartEquilibriumSolverPriorityQueue(const Partition& partition);

    /// Destroy this SmartEquilibriumSolverPriorityQueue instance.
    virtual ~SmartEquilibriumSolverPriorityQueue();

    /// Learn how to perform a full equilibrium calculation (with tracking).
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses).
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void;

    /// Output clusters created during the ODML algorithm.
    auto outputInfo() const -> void;
private:

    /// The chemical potentials at the calculated equilibrium state
    ChemicalVector u;

    /// The storage for matrix du/db, du/dT, and du/dP
    Matrix dudb, dudT, dudP;

    /// The storage for vector u(iprimary)
    Vector um;

    /// The storage for matrix Mbe = inv(u(iprimary)) * du(iprimary)/db.
    Matrix Mbe;

    /// The storage for vector MT = inv(u(iprimary)) * du(iprimary)/dT and MP = inv(u(iprimary)) * du(iprimary)/dP
    Vector MT, MP;

    struct Record
    {
        /// The temperature of the equilibrium state (in units of K).
        double T;

        /// The pressure of the equilibrium state (in units of Pa).
        double P;
        
        /// The amounts of elements in the equilibrium state (in units of mol).
        Vector be;

        /// The calculated equilibrium state at `T`, `P`, `be`.
        ChemicalState state;

        /// The chemical properties at the calculated equilibrium state.
        ChemicalProperties properties;

        /// The sensitivity derivatives at the calculated equilibrium state.
        EquilibriumSensitivity sensitivity;

        /// The matrix used to compute relative change of chemical potentials due to change in `be`.
        Matrix Mb;

        /// The vectors used to compute relative change of chemical potentials due to change in `T` and `P`.
        Vector MT, MP;

        // Indices of the major species in the chemical state. 
        VectorXi imajor;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<Record> database;

    // The queue with priority indices, in which equilibrium states of the `database` must be considered
    std::deque<Index> priority;

    // The queue with the rank of each equilibrium states in the `database`
    std::deque<Index> ranking;
};

} // namespace Reaktoro
