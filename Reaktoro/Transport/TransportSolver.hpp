// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

namespace Reaktoro {

// Forward declarations (classes)
class Mesh;

// Forward declarations (structs)
struct TransportResult;
struct TransportOptions;

/// A class for solving advection-diffusion equations.
/// Eq: du/dt + d(v*u)/dx - d/dx(D*du/dx) = q
///     u - concentration
///     v - velocity
///     D - diffusion coefficient
///     q - source rate
class TransportSolver
{
public:
    /// Construct a default TransportSolver instance.
    TransportSolver();

    /// Construct a copy of a TransportSolver instance.
    TransportSolver(const TransportSolver& other);

    /// Destroy this TransportSolver instance.
    virtual ~TransportSolver();

    /// Assign a copy of an TransportSolver instance
    auto operator=(TransportSolver other) -> TransportSolver&;

    /// Set the mesh for the numerical solution of the transport problem.
    auto setMesh(const Mesh& mesh) -> void;

    /// Set the velocity for the transport problem.
    /// @param val The velocity (in m/s)
    auto setVelocity(double val) -> void;

    /// Set the diffusion coefficient for the transport problem.
    /// @param val The diffusion coefficient (in m^2/s)
    auto setDiffusionCoeff(double val) -> void;

    /// Set the value of the variable on the boundary.
    /// @param val The boundary value for the variable (same unit considered for u).
    auto setBoundaryValue(double val) -> void;

    /// Set the time step for the numerical solution of the transport problem.
    auto setTimeStep(double val) -> void;

    /// Set options for the transport calculation.
    auto setOptions(const TransportOptions& options) -> void;

    /// Initialize the transport solver before method @ref step is executed.
    /// Setup coefficient matrix of the diffusion problem and factorize.
    auto initialize() -> void;

    /// Perform one transport time step calculation.
    /// This method performs one time step in the solution of the transport
    /// equation using an explicit scheme in time for advection and implicit
    /// scheme for diffusion.
    /// @param[in,out] u The solution vector
    /// @param q The source rates vector ([same unit considered for u]/m)
    auto step(VectorRef u, VectorConstRef q) -> TransportResult;

    /// Perform one transport time step calculation.
    /// @param[in,out] u The solution vector
    auto step(VectorRef u) -> TransportResult;

    /// Return the result of the last transport time step calculation.
    auto result() const -> const TransportResult&;

    /// Return the mesh.
    auto mesh() const -> const Mesh&;

    /// Return the last time step length used.
    auto timeStep() const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
