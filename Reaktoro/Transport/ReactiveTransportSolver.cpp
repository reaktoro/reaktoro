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

namespace Reaktoro {

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem& system)
: system_(system), equilibrium_solver(system), smart_equilibrium_solver(system)
{
    setBoundaryState(ChemicalState(system));
}

auto ReactiveTransportSolver::setOptions(const ReactiveTransportOptions& options) -> void
{
    equilibrium_solver.setOptions(options.equilibrium);
    smart_equilibrium_solver.setOptions(options.smart_equilibrium);
}

auto ReactiveTransportSolver::setMesh(const Mesh& mesh) -> void
{
    transportsolver.setMesh(mesh);
}

auto ReactiveTransportSolver::setVelocity(double val) -> void
{
    transportsolver.setVelocity(val);
}

auto ReactiveTransportSolver::setDiffusionCoeff(double val) -> void
{
    transportsolver.setDiffusionCoeff(val);
}

auto ReactiveTransportSolver::setBoundaryState(const ChemicalState& state) -> void
{
    bbc = state.elementAmounts();
}

auto ReactiveTransportSolver::setTimeStep(double val) -> void
{
    transportsolver.setTimeStep(val);
}

auto ReactiveTransportSolver::system() const -> const ChemicalSystem&
{
    system_;
}

auto ReactiveTransportSolver::output() -> ChemicalOutput
{
    outputs.emplace_back(ChemicalOutput(system_));
    return outputs.back();
}

auto ReactiveTransportSolver::initialize() -> void
{
    // Initialize mesh and corresponding amount e
    const Mesh& mesh = transportsolver.mesh();
    const Index num_elements = system_.numElements();
    const Index num_cells = mesh.numCells();

    // Initialize amount of elements in fluid and solid phases
    bf.resize(num_cells, num_elements);
    bs.resize(num_cells, num_elements);
    b.resize(num_cells, num_elements);

    // Initialize equilibrium solver based on the parameter
    transportsolver.setOptions();
    transportsolver.initialize();
}

auto ReactiveTransportSolver::step(ChemicalField& field) -> ReactiveTransportResult
{
    // The result of the reactive transport step
    ReactiveTransportResult rt_result;

    // Auxiliary variables
    const auto& mesh = transportsolver.mesh();
    const auto& num_elements = system_.numElements();
    const auto& num_cells = mesh.numCells();
    const auto& ifs = system_.indicesFluidSpecies();
    const auto& iss = system_.indicesSolidSpecies();

    // Open the the file for outputting chemical states
    for(auto output : outputs)
    {
        output.suffix("-" + std::to_string(steps));
        output.open();
    }

    //---------------------------------------------------------------------------
    // Step 1: Perform a time step transport calculation for each fluid element
    //---------------------------------------------------------------------------
    tic();

    // Collect the amounts of elements in the solid and fluid species
    for(Index icell = 0; icell < num_cells; ++icell)
    {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Left boundary condition cell
    Index icell_bc = 0;
    const auto phi_bc = field[icell_bc].properties().fluidVolume().val;

    // Transport the elements in the fluid species
    for(Index ielement = 0; ielement < num_elements; ++ielement)
    {
        // Scale BC with a porosity of the boundary cell
        transportsolver.setBoundaryValue(phi_bc * bbc[ielement]);
        transportsolver.step(bf.col(ielement));
    }

    // Sum the amounts of elements distributed among fluid and solid species
    b.noalias() = bf + bs;

    toc() >> profiling.time_transport;

    //---------------------------------------------------------------------------
    // Step 1: Perform a time step transport calculation for each fluid element
    //---------------------------------------------------------------------------
    tic();

    if(options.use_smart_equilibrium_solver)
    {
        SmartEquilibriumResult eqres;

        // Ensure the result and profiling info of each cell's smart equilibrium calculation can be saved.
        ifprofiling( profiling.smart_equilibrium_result_at_cell.resize(num_cells); );
        ifprofiling( profiling.smart_equilibrium_profiling_at_cell.resize(num_cells); );

        for(Index icell = 0; icell < num_cells; ++icell)
        {
            const auto T = field[icell].temperature();
            const auto P = field[icell].pressure();

            // Solve with a smart equilibrium solver
            eqres = smart_equilibrium_solver.solve(field[icell], T, P, b.row(icell));

            // Store the result and profiling info of this cell's smart equilibrium calculation.
            ifprofiling( profiling.smart_equilibrium_result_at_cell[icell] = eqres; );
            ifprofiling( profiling.smart_equilibrium_profiling_at_cell[icell] = smart_equilibrium_solver.profiling(); );
        }
    }
    else
    {
        EquilibriumResult eqres;

        // Ensure the result and profiling info of each cell's equilibrium calculation can be saved.
        ifprofiling( profiling.equilibrium_result_at_cell.resize(num_cells); );
        ifprofiling( profiling.equilibrium_profiling_at_cell.resize(num_cells); );

        for(Index icell = 0; icell < num_cells; ++icell)
        {
            const auto T = field[icell].temperature();
            const auto P = field[icell].pressure();

            // Solve with a conventional equilibrium solver
            eqres = equilibrium_solver.solve(field[icell], T, P, b.row(icell));

            // Store the result and profiling info of this cell's smart equilibrium calculation.
            ifprofiling( profiling.equilibrium_result_at_cell[icell] = eqres; );
            ifprofiling( profiling.equilibrium_profiling_at_cell[icell] = equilibrium_solver.profiling(); );
        }
    }

    toc() >> profiling.time_equilibrium;


    // for(Index icell = 0; icell < num_cells; ++icell)
    // {
    //     const double T = field[icell].temperature();
    //     const double P = field[icell].pressure();

    //     // Start profiling equlibrium
    //     ifprofiling( start = time(); );

    //     if(options.use_smart_equilibrium_solver)
    //     {
    //         // Solve with a smart equilibrium solver
    //         rt_result.equilibrium += smart_equilibrium_solver.solve(field[icell], T, P, b.row(icell));

    //         // End profiling for the equilibrium calculations (accumulate cell-wise)
    //         ifprofiling( rt_result.eq_time += elapsed(start); );

    //         // Update the time spend for either for learning or estimating
    //         if(rt_result.equilibrium.smart.succeeded)
    //             rt_result.equilibrium_smart_successfull_cells.push_back(icell);
    //     }
    //     else
    //     {
    //         // Solve with a conventional equilibrium solver
    //         rt_result.equilibrium += equilibrium_solver.solve(field[icell], T, P, b.row(icell));

    //         // End profiling for the conventional equilibrium calculations (accumulate cell-wise)
    //         ifprofiling( rt_result.eq_time += elapsed(start); );
    //     }
    // }

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

} // namespace Reaktoro
