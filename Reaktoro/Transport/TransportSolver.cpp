// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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


#include "TransportSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Transport/Mesh.hpp>
#include <Reaktoro/Transport/TransportOptions.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>
#include <Reaktoro/Transport/TridiagonalMatrix.hpp>

namespace Reaktoro {

struct TransportSolver::Impl
{
    /// The mesh describing the discretization of the domain.
    Mesh mesh;

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
    ArrayXr phi;

    /// The previous state of the variables.
    ArrayXr u0;

    /// The options for the transport calculations.
    TransportOptions options;

    /// The result of the last transport time step calculation.
    TransportResult result;

    /// Construct a default TransportSolver::Impl instance.
    Impl()
    {
    }

    /// Set the mesh for the numerical solution of the transport problem.
    auto setMesh(const Mesh& mesh) -> void
    {
        this->mesh = mesh;
    }

    /// Set the velocity for the transport problem.
    auto setVelocity(double val) -> void
    {
        velocity = val;
    }

    /// Set the diffusion coefficient for the transport problem.
    auto setDiffusionCoeff(double val) -> void
    {
        diffusion = val;
    }

    /// Set the value of the variable on the boundary.
    auto setBoundaryValue(double val) -> void
    {
        ul = val;
    }

    /// Set the time step for the numerical solution of the transport problem.
    auto setTimeStep(double val) -> void
    {
        dt = val;
    }

    /// Set the options for the transport solver.
    auto setOptions(const TransportOptions& options) -> void
    {
        this->options = options;
    }

    /// Initialize the transport solver before method @ref step is executed.
    auto initialize() -> void
    {
        const auto dx = mesh.dx();
        const auto alpha = diffusion*dt/(dx * dx);
        const auto beta = velocity*dt/dx;
        const auto num_cells = mesh.numCells();
        const auto icell0 = 0;
        const auto icelln = num_cells - 1;

        tic(ASSEMBLY_STEP);

        // Initialize A vector
        A.resize(num_cells);

        // Initialize coefficients of the system's matrix
        double a(0.0), b(0.0), c(0.0);

        // Assemble the coefficient matrix A for the interior cells
        for(Index icell = 1; icell < icelln; ++icell)
        {
            // Depending on the initialized FV scheme, initialize coefficients
            switch(options.finite_volume_method)
            {
                case FiniteVolumeMethod::FullImplicit:
                    a = -beta - alpha;
                    b = 1 + beta + 2 * alpha;
                    c = -alpha;
                    break;
            }
            A.row(icell) << a, b, c;
        }

        // Initialize cells on the left and right BCs
        switch(options.finite_volume_method)
        {
            case FiniteVolumeMethod::FullImplicit:
                // Assemble the coefficient matrix A for the boundary cells
                A.row(icell0) << 0.0, 1.0 + alpha + beta, -alpha;   // left boundary (flux = v * ul)
                A.row(icelln) << - beta, 1.0 + beta, 0.0;           // right boundary (free)

                break;
        }

        // Factorize A into LU factors for future uses in method step
        A.factorize();

        result.timing.matrix_equation_assembly = toc(ASSEMBLY_STEP);
    }

    /// Perform one transport time step calculation.
    auto step(ArrayXrRef u, ArrayXrConstRef q) -> TransportResult
    {
        tic(TRANSPORT_STEP);

        // Reset the result of the last transport calculation
        result = {};

        // Solving advection problem with time explicit approach
        const auto dx = mesh.dx();
        const auto beta = velocity*dt/dx;
        const auto icell0 = 0;

        switch(options.finite_volume_method)
        {
            case FiniteVolumeMethod::FullImplicit:
                // Handle the left boundary cell
                u[icell0] += beta * ul; // left boundary (flux = v * ul)
                break;
        }

        // Add the source contribution
        u += dt * q;

        // Solving the diffusion problem with time implicit approach
        timeit( A.solve(u), result.timing.matrix_equation_solve= );

        result.timing.step = toc(TRANSPORT_STEP);

        return result;
    }

    /// Perform one transport time step calculation.
    auto step(ArrayXrRef u) -> TransportResult
    {
        return step(u, ArrayXr::Zero(u.size()));
    }
};

TransportSolver::TransportSolver()
: pimpl(new Impl())
{}

TransportSolver::TransportSolver(const TransportSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

TransportSolver::~TransportSolver()
{}

auto TransportSolver::operator=(TransportSolver other) -> TransportSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto TransportSolver::setMesh(const Mesh& mesh) -> void
{
    pimpl->setMesh(mesh);
}

auto TransportSolver::setVelocity(double val) -> void
{
    pimpl->setVelocity(val);
}

auto TransportSolver::setDiffusionCoeff(double val) -> void
{
    pimpl->setDiffusionCoeff(val);
}

auto TransportSolver::setBoundaryValue(double val) -> void
{
    pimpl->setBoundaryValue(val);
}

auto TransportSolver::setTimeStep(double val) -> void
{
    pimpl->setTimeStep(val);
}

auto TransportSolver::setOptions(const TransportOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto TransportSolver::initialize() -> void
{
    pimpl->initialize();
}

auto TransportSolver::step(ArrayXrRef u, ArrayXrConstRef q) -> TransportResult
{
    return pimpl->step(u, q);
}

auto TransportSolver::step(ArrayXrRef u) -> TransportResult
{
    return pimpl->step(u);
}

auto TransportSolver::result() const -> const TransportResult&
{
    return pimpl->result;
}

auto TransportSolver::mesh() const -> const Mesh&
{
    return pimpl->mesh;
}

auto TransportSolver::timeStep() const -> double
{
    return pimpl->dt;
}

} // namespace Reaktoro
