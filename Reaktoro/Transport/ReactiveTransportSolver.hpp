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
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Transport/ReactiveTransportOptions.hpp>
#include <Reaktoro/Transport/ReactiveTransportResult.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalField;
class ChemicalOutput;
class ChemicalSystem;
class ChemicalState;
class Mesh;
class Partition;

/// Use this class for solving reactive transport problems.
class ReactiveTransportSolver
{
public:
    /// Construct a ReactiveTransportSolver instance with given chemical system.
    explicit ReactiveTransportSolver(const ChemicalSystem& system);

    /// Construct a ReactiveTransportSolver instance with given partition.
    explicit ReactiveTransportSolver(const Partition& partition);

    /// Construct a copy of a ReactiveTransportSolver instance.
    ReactiveTransportSolver(const ReactiveTransportSolver& other);

    /// Destroy this ReactiveTransportSolver instance.
    virtual ~ReactiveTransportSolver();

    /// Assign a copy of an ReactiveTransportSolver instance
    auto operator=(ReactiveTransportSolver other) -> ReactiveTransportSolver&;

    /// Set the options for the reactive transport calculations.
    auto setOptions(const ReactiveTransportOptions& options) -> void;

    /// Initialize the mesh discretizing the computational domain for reactive transport.
    auto setMesh(const Mesh& mesh) -> void;

    /// Initialize the velocity of the reactive transport model.
    auto setVelocity(double val) -> void;

    /// Initialize the diffusion of the reactive transport model.
    auto setDiffusionCoeff(double val) -> void;

    /// Initialize boundary conditions of the reactive transport model.
    /// Method initializes the boundary condition of the reactive transport model
    /// by setting there equilibrated ChemicalState
    auto setBoundaryState(const ChemicalState& state) -> void;

    /// Initialize time step of the reactive transport sequential algorithm.
    auto setTimeStep(double val) -> void;

    /// Return a ChemicalOutput instance.
    /// The returned ChemicalOutput instance must be properly configured
    /// before the method ReactiveTransportSolver::step is called.
    /// Changes in this ChemicalOutput instance are observed by this
    /// ReactiveTransportSolver object.
    auto output() -> ChemicalOutput;

    /// Initialize the reactive transport solver.
    /// This method should be called before the first call to method ReactiveTransportSolver::step.
    auto initialize() -> void;

    /// Perform one time step of a reactive transport calculation.
    auto step(ChemicalField& field) -> ReactiveTransportResult;

    /// Return the result of the last reactive transport time step calculation.
    auto result() const -> const ReactiveTransportResult&;

    /// Return the chemical system used in the reactive transport solver.
    auto system() const -> const ChemicalSystem&;

    /// Return the last time step length used.
    auto timeStep() const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
