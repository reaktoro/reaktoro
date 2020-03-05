// This file is part of Reaktoro (https://reaktoro.org).
//
// Reaktoro is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// Reaktoro is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ReactiveTransportSolver.hpp"

// C++ includes
#include <algorithm>
#include <fstream>
#include <iomanip>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Transport/Mesh.hpp>

namespace Reaktoro {

struct ReactiveTransportSolver::Impl
{
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem system;

    /// The solver for solving the transport equations
    TransportSolver transport_solver;

    /// The options for the reactive transport calculations.
    ReactiveTransportOptions options;

    /// The result information of the last reactive transport time step calculation.
    ReactiveTransportResult result;

    /// The equilibrium solver using conventional Gibbs energy minimization approach.
    EquilibriumSolver equilibrium_solver;

    /// The equilibrium solver using a smart on-demand learning strategy.
    SmartEquilibriumSolver smart_equilibrium_solver;

    /// The list of chemical output objects
    std::vector<ChemicalOutput> outputs;

    /// The amounts of fluid elements on the boundary.
    Vector bbc;

    /// The amounts of a fluid element on each cell of the mesh.
    Matrix bf;

    /// The amounts of a solid element on each cell of the mesh.
    Matrix bs;

    /// The amounts of an element on each cell of the mesh.
    Matrix b;

    /// The current number of steps in the solution of the reactive transport equations.
    Index steps = 0;

    /// Name of the file and folder with a status output
    std::string folder;

    /// Construct a ReactiveTransportSolver::Impl instance.
    Impl(const ChemicalSystem& system)
    : system(system), equilibrium_solver(system), smart_equilibrium_solver(system)
    {
        setBoundaryState(ChemicalState(system));
    }

    /// Set the options for the reactive transport calculations.
    auto setOptions(const ReactiveTransportOptions& options) -> void
    {
        this->options = options;
        equilibrium_solver.setOptions(options.equilibrium);
        smart_equilibrium_solver.setOptions(options.smart_equilibrium);
    }

    /// Initialize the mesh discretizing the computational domain for reactive transport.
    auto setMesh(const Mesh& mesh) -> void
    {
        transport_solver.setMesh(mesh);
    }

    /// Initialize the velocity of the reactive transport model.
    auto setVelocity(double val) -> void
    {
        transport_solver.setVelocity(val);
    }

    /// Initialize the diffusion of the reactive transport model.
    auto setDiffusionCoeff(double val) -> void
    {
        transport_solver.setDiffusionCoeff(val);
    }

    /// Initialize boundary conditions of the reactive transport model.
    auto setBoundaryState(const ChemicalState& state) -> void
    {
        bbc = state.elementAmounts();
    }

    /// Initialize time step of the reactive transport sequential algorithm.
    auto setTimeStep(double val) -> void
    {
        transport_solver.setTimeStep(val);
    }

    /// Add the output to the reactive transport modelling.
    auto output() -> ChemicalOutput
    {
        outputs.emplace_back(ChemicalOutput(system));
        return outputs.back();
    }

    /// Initialize the reactive transport solver.
    auto initialize() -> void
    {
        // Initialize mesh and corresponding amount e
        const Mesh& mesh = transport_solver.mesh();
        const Index num_elements = system.numElements();
        const Index num_cells = mesh.numCells();

        // Initialize amount of elements in fluid and solid phases
        bf.resize(num_cells, num_elements);
        bs.resize(num_cells, num_elements);
        b.resize(num_cells, num_elements);

        // Initialize equilibrium solver based on the parameter
        transport_solver.setOptions(options.transport);
        transport_solver.initialize();
    }

    /// Perform one time step of a reactive transport calculation.
    auto step(ChemicalField& field) -> ReactiveTransportResult
    {
        // The result of the reactive transport step
        ReactiveTransportResult rt_result;

        // Auxiliary variables
        const auto& mesh = transport_solver.mesh();
        const auto& num_elements = system.numElements();
        const auto& num_cells = mesh.numCells();
        const auto& ifs = system.indicesFluidSpecies();
        const auto& iss = system.indicesSolidSpecies();

        // Open the the file for outputting chemical states
        for(auto output : outputs)
        {
            output.suffix("-" + std::to_string(steps));
            output.open();
        }

        //---------------------------------------------------------------------------
        // Step 1: Perform a time step transport calculation for each fluid element
        //---------------------------------------------------------------------------
        tic(TRANSPORT_STEP);

        // Collect the amounts of elements in the solid and fluid species
        for(Index icell = 0; icell < num_cells; ++icell)
        {
            bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
            bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
        }

        // Left boundary condition cell
        Index icell_bc = 0;
        // Get porosity of the left boundary cell
        const auto phi_bc = field[icell_bc].properties().fluidVolume().val / field[icell_bc].properties().volume().val;

        // Ensure the result of each fluid element transport calculation can be saved
        result.transport_of_element.resize(num_elements);

        // Transport the elements in the fluid species
        for(Index ielement = 0; ielement < num_elements; ++ielement)
        {
            // Scale BC with a porosity of the boundary cell
            transport_solver.setBoundaryValue(phi_bc * bbc[ielement]);
            transport_solver.step(bf.col(ielement));

            // Save the result of this element transport calculation.
            result.transport_of_element[ielement] = transport_solver.result();
        }

        // Sum the amounts of elements distributed among fluid and solid species
        b.noalias() = bf + bs;

        result.timing.transport = toc(TRANSPORT_STEP);

        //---------------------------------------------------------------------------
        // Step 2: Perform a time step equilibrium calculation for each cell
        //---------------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP);

        if(options.use_smart_equilibrium_solver)
        {
            // Ensure the result of each cell's smart equilibrium calculation can be saved
            result.smart_equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = field[icell].temperature();
                const auto P = field[icell].pressure();

                // Solve with a smart equilibrium solver
                smart_equilibrium_solver.solve(field[icell], T, P, b.row(icell));

                // Save the result of this cell's smart equilibrium calculation
                result.smart_equilibrium_at_cell[icell] = smart_equilibrium_solver.result();
            }
        }
        else
        {
            // Ensure the result of each cell's equilibrium calculation can be saved.
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = field[icell].temperature();
                const auto P = field[icell].pressure();

                // Solve with a conventional equilibrium solver
                equilibrium_solver.solve(field[icell], T, P, b.row(icell));

                // Save the result of this cell's smart equilibrium calculation.
                result.equilibrium_at_cell[icell] = equilibrium_solver.result();
            }
        }

        result.timing.equilibrium = toc(EQUILIBRIUM_STEP);

        // Update the output files with the chemical state of every cell
        for(Index icell = 0; icell < num_cells; ++icell)
            for(auto output : outputs)
                output.update(field[icell], icell);

        // Output chemical states in the output files
        for(auto output : outputs)
            output.close();

        // Increment the current number of reactive transport steps
        ++steps;

        return rt_result;
    }
};

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{
}

ReactiveTransportSolver::ReactiveTransportSolver(const ReactiveTransportSolver& other)
: pimpl(new Impl(*other.pimpl))
{
}

ReactiveTransportSolver::~ReactiveTransportSolver()
{
}

auto ReactiveTransportSolver::operator=(ReactiveTransportSolver other) -> ReactiveTransportSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ReactiveTransportSolver::setOptions(const ReactiveTransportOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto ReactiveTransportSolver::setMesh(const Mesh& mesh) -> void
{
    pimpl->setMesh(mesh);
}

auto ReactiveTransportSolver::setVelocity(double val) -> void
{
    pimpl->setVelocity(val);
}

auto ReactiveTransportSolver::setDiffusionCoeff(double val) -> void
{
    pimpl->setDiffusionCoeff(val);
}

auto ReactiveTransportSolver::setBoundaryState(const ChemicalState& state) -> void
{
    pimpl->setBoundaryState(state);
}

auto ReactiveTransportSolver::setTimeStep(double val) -> void
{
    pimpl->setTimeStep(val);
}

auto ReactiveTransportSolver::output() -> ChemicalOutput
{
    return pimpl->output();
}

auto ReactiveTransportSolver::initialize() -> void
{
    pimpl->initialize();
}

auto ReactiveTransportSolver::step(ChemicalField& field) -> ReactiveTransportResult
{
    return pimpl->step(field);
}

auto ReactiveTransportSolver::result() const -> const ReactiveTransportResult&
{
    return pimpl->result;
}

auto ReactiveTransportSolver::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ReactiveTransportSolver::timeStep() const -> double
{
    return pimpl->transport_solver.timeStep();
}

} // namespace Reaktoro
