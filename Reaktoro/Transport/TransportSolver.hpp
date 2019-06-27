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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Transport/Mesh.hpp>

namespace Reaktoro {

//class BoundaryState
//{
//public:
//    BoundaryState(const ChemicalState& state);
//
//private:
//};
//
//class BoundaryFlux
//{
//public:
//    BoundaryFlux(const ChemicalState& state);
//
//private:
//};

enum FiniteVolumeMethod{
    FullImplicit = 1,
    ImpliciteExpilcit = 2,
    FluxLimitersImplicitExplicit = 3
};

/// A class that defines a Tridiagonal Matrix used on TransportSolver.
/// it stores data in a Eigen::VectorXd like, M = {a[0][0], a[0][1], a[0][2],
///                                                a[1][0], a[1][1], a[1][2],
///                                                a[2][0], a[2][1], a[2][2]}
class TridiagonalMatrix
{
public:
    TridiagonalMatrix() : TridiagonalMatrix(0) {}

    TridiagonalMatrix(Index size) : m_size(size), m_data(size * 3) {}

    auto size() const -> Index { return m_size; }

    auto data() -> VectorRef { return m_data; }

    auto data() const -> VectorConstRef { return m_data; }

    auto row(Index index) -> VectorRef { return m_data.segment(3 * index, 3); }

    auto row(Index index) const -> VectorConstRef { return m_data.segment(3 * index, 3); }

    auto a() -> VectorStridedRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto a() const -> VectorConstRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto b() -> VectorStridedRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto b() const -> VectorConstRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto c() -> VectorStridedRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto c() const -> VectorConstRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto resize(Index size) -> void;

    /// Factorize the tridiagonal matrix to help with linear system A x = b.
    auto factorize() -> void;

    /// Solve a linear system Ax = b with LU decomposition.
    auto solve(VectorRef x, VectorConstRef b) const -> void;

    /// Solve a linear system Ax = b with LU decomposition, using x as the unknown and it's.
    /// old values as the vector b.
    auto solve(VectorRef x) const -> void;

    operator Matrix() const;

private:
    /// The size of the tridiagonal matrix
    Index m_size;

    /// The coefficients
    Vector m_data;
};

/// A class for solving advection-diffusion problem.
/// Eq: du/dt + v*du/dx = D*d2u/dx2
///     u - amount
///     v - velocity
///     D - diffusion coefficient
class TransportSolver
{
public:
    /// Construct a default TransportSolver instance.
    TransportSolver();

    /// Set the mesh for the numerical solution of the transport problem.
    auto setMesh(const Mesh& mesh) -> void { mesh_ = mesh; }

    /// Set the velocity for the transport problem.
    /// @param val The velocity (in m/s)
    auto setVelocity(double val) -> void { velocity = val; }

    /// Set the diffusion coefficient for the transport problem.
    /// @param val The diffusion coefficient (in m^2/s)
    auto setDiffusionCoeff(double val) -> void { diffusion = val; }

    /// Set the value of the variable on the boundary.
    /// @param val The boundary value for the variable (same unit considered for u).
    auto setBoundaryValue(double val) -> void { ul = val; };

    /// Set the time step for the numerical solution of the transport problem.
    auto setTimeStep(double val) -> void { dt = val; }

    /// Set the time step for the numerical solution of the transport problem.
    auto setOptions() -> void { options.fvm = FullImplicit; }

    /// Return the mesh.
    auto mesh() const -> const Mesh& { return mesh_; }

    /// Initialize the transport solver before method @ref step is executed.
    /// Setup coefficient matrix of the diffusion problem and factorize.
    auto initialize() -> void;

    /// Step the transport solver.
    /// This method solve one step of the transport solver equation, using an explicit approach for
    /// advection and total implicit for diffusion. The amount resulted from the advection it is
    /// passed to diffusion problem as a "source".
    /// @param[in,out] u The solution vector
    /// @param q The source rates vector ([same unit considered for u]/m)
    auto step(VectorRef u, VectorConstRef q) -> void;

    /// Step the transport solver.
    /// @param[in,out] u The solution vector
    auto step(VectorRef u) -> void;

private:
    /// The mesh describing the discretization of the domain.
    Mesh mesh_;

    /// The time step used to solve the transport problem (in s).
    double dt = 0.0;

    /// The velocity in the transport problem (in m/s).
    double velocity = 0.0;

    /// The diffusion coefficient in the transport problem (in m^2/s).
    double diffusion = 0.0;

    /// The value of the variable on the left boundary.
    double ul;

    /// The coefficient matrix from the discretized transport equation
    TridiagonalMatrix A;

    /// The flux limiters at each cell.
    Vector phi;

    /// The previous state of the variables.
    Vector u0;

    struct Options{

        /// Flag for the FV scheme
        FiniteVolumeMethod fvm = FullImplicit;
    };

    Options options;
};

} // namespace Reaktoro
