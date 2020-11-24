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
#include <Reaktoro/ODML/ClusterConnectivity.hpp>
#include <Reaktoro/ODML/PriorityQueue.hpp>

namespace Reaktoro {
    
/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolverNN : public SmartEquilibriumSolverBase
{
public:

    /// Construct an SmartEquilibriumSolverNN instance with given chemical system.
    SmartEquilibriumSolverNN(const ChemicalSystem& system);

    /// Construct an SmartEquilibriumSolverNN instance with given partition of the chemical system.
    SmartEquilibriumSolverNN(const Partition& partition);

    /// Destroy this SmartEquilibriumSolverNN instance.
    virtual ~SmartEquilibriumSolverNN();

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void;

    /// Output clusters created during the ODML algorithm
    auto outputInfo() const -> void;
private:

    // The structure to store the information about learned equilibrium
    struct TreeNode
    {
        Vector be;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
    };
    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<TreeNode> tree;
};

} // namespace Reaktoro
