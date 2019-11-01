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
    Vector phi;

    /// The previous state of the variables.
    Vector u0;

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

        tic(0);

        // Initialize A vector
        A.resize(num_cells);

        // If the Flux limiters scheme is considered, initialize phi vector
        if(options.finite_volume_method == FiniteVolumeMethod::FluxLimitersImplicitExplicit)
        {
            Assert(velocity * dt / dx < 0.5,
                "Could not run reactive-transport calculation reliably.",
                "The CFL number ( = velocity * dt / dx ) must be less then 0.5. "
                "Try to decrease the time step or coarsen the spatial discretization"
                "(increase the number of cells)");
            phi.resize(num_cells);
        }

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

                case FiniteVolumeMethod::FluxLimitersImplicitExplicit:
                    a = -alpha;
                    b = 1 + 2*alpha;
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

            case FiniteVolumeMethod::FluxLimitersImplicitExplicit:

                // Assemble the coefficient matrix A for the boundary cells
                // Our derivation
                // A.row(icell0) << 0.0, 1.0 + alpha, -alpha;            // forward difference approximation with first order error
                // A.row(icelln) << 0.0, 1.0, 0.0;                       // right boundary (free)

                // Allan version ealier
                // A.row(icell0) << 0.0, 1.0 + alpha, -alpha;            // forward difference approximation with second order error
                // A.row(icelln) << -alpha, 1.0 + alpha, 0.0;            // d/dx = 0 (zero flux) at the right boundary (to be fixed)

                // ESSS
                A.row(icell0) << 0.0, 1.0 + 4.5*alpha, -1.5*alpha;    // forward difference approximation with second order error
                A.row(icelln) << -alpha, 1.0 + alpha, 0.0;            // d/dx = 0 (zero flux) at the right boundary

                break;
        }

        // Factorize A into LU factors for future uses in method step
        A.factorize();

        toc(0, result.timing.matrix_equation_assembly);

    }

    /// Perform one transport time step calculation.
    auto step(VectorRef u, VectorConstRef q) -> TransportResult
    {
        tic(0);

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

            case FiniteVolumeMethod::FluxLimitersImplicitExplicit:

                const auto num_cells = mesh.numCells();
                const auto icelln = num_cells - 1;

                u0 = u;

                phi[0] = 2.0; //  this is very important to ensure correct flux limiting behavior for boundary cell.

                // Calculate the flux limiters in the interior cells
                for(Index icell = 1; icell < icelln; ++icell)
                {
                    // Calculate the variation index `r = (uP - uW)/(uE - uP)` on current cell
                    const double r = (u[icell] - u[icell - 1])/(u[icell + 1] - u[icell]);

                    // Calculate the flux limiter phi based on the superbee limiter (https://en.wikipedia.org/wiki/Flux_limiter)
                    phi[icell] = std::max(0.0, std::max(std::min(2 * r, 1.0), std::min(r, 2.0)));
                }

                // Compute advection contributions to u for the interior cells
                for(Index icell = 1; icell < icelln; ++icell)
                {
                    const double phiW = phi[icell - 1];
                    const double phiP = phi[icell];
                    const double aux = 1.0 + 0.5 * (phiP - phiW);

                    const double uW = u0[icell - 1];
                    const double uP = u0[icell];
                    u[icell] += aux*beta * (uW - uP);
                }

                // Handle the left boundary cell
                const double aux = 1 + 0.5 * phi[0];
                u[icell0] += aux * beta * (ul - u0[0]) + (3.0*diffusion*ul*dt/(dx*dx)); // prescribed amount on the wall and approximation derived by forward difference approximation with second order error
                // Handle the right boundary cell
                u[icelln] += beta * (u0[icelln - 1] - u0[icelln]); // du/dx = 0 at the right boundary
                break;
        }

        // Add the source contribution
        u += dt * q;

        // Solving the diffusion problem with time implicit approach
        timeit( A.solve(u), result.timing.matrix_equation_solve= );

        toc(0, result.timing.step);

        return result;
    }

    /// Perform one transport time step calculation.
    auto step(VectorRef u) -> TransportResult
    {
        return step(u, zeros(u.size()));
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

auto TransportSolver::step(VectorRef u, VectorConstRef q) -> TransportResult
{
    return pimpl->step(u, q);
}

auto TransportSolver::step(VectorRef u) -> TransportResult
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

//auto TransportSolver::step(VectorRef u, VectorConstRef q) -> void
//{
//    // TODO: Implement Kurganov-Tadmor method as detailed in their 2000 paper (not as in Wikipedia)
//    const double dx = mmesh.dx();
//    const double num_cells = mmesh.numCells();
//    const double alpha = velocity*dt/dx;
//    const double beta = diffusion*dt/(dx * dx);
//
//    u0 = u;
//    A.resize(num_cells);
//    phi.resize(num_cells);
//
//    // Calculate the flux limiters in the interior cells
//    for(Index i = 1; i < num_cells - 1; ++i)
//    {
//        // Calculate the variation index `r = (uP - uW)/(uE - uP)` on current cell
//        const double r = (u[i] - u[i - 1])/(u[i + 1] - u[i]);
//
//        // Calculate the flux limiter phi based on the superbee limiter (https://en.wikipedia.org/wiki/Flux_limiter)
//        phi[i] = std::max(0.0, std::max(std::min(2 * r, 1.0), std::min(r, 2.0)));
//    }
//
//    // Assemble the coefficient matrix A for the interior cells
//    for(Index icell = 1; icell < num_cells - 1; ++icell)
//    {
//        const double phiW = phi[icell - 1];
//        const double phiP = phi[icell];
//        const double aux = 1.0 + 0.5 * (phiP - phiW);
//        const double a = -beta;
//        const double b = 1 + 2*beta;
//        const double c = -beta;
//        A.row(icell) << a, b, c;
//
//        const double uW = u0[icell - 1];
//        const double uP = u0[icell];
//        u[icell] += aux*alpha * (uW - uP);
//    }
//
//    // Assemble the coefficient matrix A for the boundary cells
//    A.row(0) << 0.0, 1.0 + beta, -beta;
//    A.row(num_cells - 1) << -beta, 1.0 + beta, 0.0;
//
//    u[0] += alpha * (ul - u0[0]);
//    u[num_cells - 1] += alpha * (u0[num_cells - 2] - u0[num_cells - 1]);
//
//    u += dt * q;
//
//    A.factorize();
//    A.solve(u);
//}

} // namespace Reaktoro
