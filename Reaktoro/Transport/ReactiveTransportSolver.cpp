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

// C++ includes
#include <algorithm>
#include <fstream>
#include <iomanip>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Dense>

// Reaktoro includes
#include "Reaktoro/Common/Exception.hpp"
#include "Reaktoro/Transport/ReactiveTransportSolver.hpp"
#include "Reaktoro/Common/TimeUtils.hpp"

namespace Reaktoro {

/// Implementation of a ReactiveTransportResult that collects information per one step of the reactive transport
///
ReactiveTransportResult::ReactiveTransportResult(const int & ncells, const int & nsteps, const bool& smart )
    : ncells(ncells), nsteps(nsteps), smart(smart){}
auto ReactiveTransportResult::resetTime() -> void
{

    if (smart) {
        // Reset the CPU times
        equilibrium.smart.estimate_stats.time_estimate = 0.0;
        equilibrium.smart.estimate_stats.time_search = 0;
        equilibrium.smart.estimate_stats.time_mat_vect_mult = 0;
        equilibrium.smart.estimate_stats.time_acceptance = 0;
        equilibrium.smart.learn_stats.time_learn = 0;
        equilibrium.smart.learn_stats.time_store = 0;
        equilibrium.smart.learn_stats.time_gibbs_min = 0;

        // Clear the indices of the learned states
        equilibrium.smart.learning_states_indx.clear();
    } else{
        equilibrium.stats.time_learn = 0;
    }
    eq_time = 0;
    rt_time = 0;

}

/// Implementation of a ReactiveTranportProfiler that collects information accumulated during the reactive transport
///
ReactiveTransportProfiler::ReactiveTransportProfiler(const std::string &folder,
        const std::string & file, const bool & smart)
    : folder(folder), file(file), smart(smart) {}

/// Process results collected per each step of reactive transport
auto ReactiveTransportProfiler::process(ReactiveTransportResult &rt_result) -> void {
    if (smart) {
        SmartEquilibriumResult smart_res = rt_result.equilibrium.smart;
        Statistics stats;

        // Initialize the results step-wise
        stats.time_estimate = smart_res.estimate_stats.time_estimate;
        stats.time_search = smart_res.estimate_stats.time_search;
        stats.time_mat_vect_mult = smart_res.estimate_stats.time_mat_vect_mult;
        stats.time_acceptance = smart_res.estimate_stats.time_acceptance;
        stats.time_learn = smart_res.learn_stats.time_learn;
        stats.time_store = smart_res.learn_stats.time_store;
        stats.time_gibbs_min = smart_res.learn_stats.time_gibbs_min;
        stats.tree_size = smart_res.tree_size;
        stats.learning_counter = (smart_res.learning_states_indx.empty()) ? 0 : smart_res.learning_states_indx.size();
        step_stats.emplace_back(stats);

        // Accumulate the step-wise results
        total_stats.time_estimate += smart_res.estimate_stats.time_estimate;
        total_stats.time_search += smart_res.estimate_stats.time_search;
        total_stats.time_mat_vect_mult += smart_res.estimate_stats.time_mat_vect_mult;
        total_stats.time_acceptance += smart_res.estimate_stats.time_acceptance;
        total_stats.time_learn += smart_res.learn_stats.time_learn;
        total_stats.time_store += smart_res.learn_stats.time_store;
        total_stats.time_gibbs_min += smart_res.learn_stats.time_gibbs_min;
        total_stats.tree_size += smart_res.tree_size;
        total_stats.learning_counter += stats.learning_counter ;
        total_stats.total_counter += rt_result.ncells;

        // Fill in the statuses array based on the
        auto it = smart_res.learning_states_indx.begin();

        for (Index i = 0; i < rt_result.ncells; i++) {
            if (it != smart_res.learning_states_indx.end() && *it == i){
                statuses.emplace_back(false);
                it++;
            }
            else
                statuses.emplace_back(true);
        }

    }
    else{
        EquilibriumResult res = rt_result.equilibrium;
        Statistics stats;

        // Initialize the results step-wise
        stats.time_learn = res.stats.time_learn;
        step_stats.emplace_back(stats);

        // Accumulate the step-wise results
        total_stats.time_learn += res.stats.time_learn;
    }
    // Add reactive transport and equilibrium time
    rt_times.emplace_back(rt_result.rt_time);
    eq_times.emplace_back(rt_result.eq_time);
    // Null the CUP time
    rt_result.resetTime();

    };

/// Output results collected at each step of reactive transport
auto ReactiveTransportProfiler::output(const Index &step) -> void {
    // The output stream of the data file
    std::ofstream datafile;

    // The floating-point precision in the output.
    int precision = 6;

    // Statuses for opening of the file
    auto opt = std::ofstream::out;

    // Depending on the step, open new file or append to existing one
    if (step != 0 && !file.empty()) opt |= std::ofstream::app;
    else if (step == 0 && !file.empty()) opt |= std::ofstream::trunc;

    if (smart) {
        // Document statuses
        // ----------------------------------------------------------------------
        datafile.open(folder + "/statuses.txt", opt);
        // Output statuses such that steps' go vertically and cells horizontally
        if (datafile.is_open()) {
            // Output statuses collected while stepping with ReactiveTransportSolver
            for (bool st : statuses) datafile << std::to_string(st) << "\t";
            datafile << "\n";
        }
        statuses.clear();

        // Close the data file
        datafile.close();
    }

}

/// Summarize the profiling results of reactive transport
auto ReactiveTransportProfiler::summarize(Results& results) -> void
{
    if (smart) {
        std::cout << std::setprecision(6);

        std::cout << "\nlearning time                  : " << total_stats.time_learn << " (" << total_stats.time_learn / (total_stats.time_learn + total_stats.time_estimate) * 100 << "% of total time)\n";
        std::cout << "   - gibbs minimization time   : " << total_stats.time_gibbs_min << " (" << total_stats.time_gibbs_min / total_stats.time_learn * 100 << "% of learning time)\n";
        std::cout << "   - store time                : " << total_stats.time_store << " (" << total_stats.time_store / total_stats.time_learn * 100 << "% of learning time)\n";
        std::cout << "estimating time                : " << total_stats.time_estimate  << " (" << total_stats.time_estimate / (total_stats.time_learn + total_stats.time_estimate) * 100 << "% of total time)\n";
        std::cout << "   - search time               : " << total_stats.time_search << " (" << total_stats.time_search / total_stats.time_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - matrix-vector mult. time  : " << total_stats.time_mat_vect_mult << " (" << total_stats.time_mat_vect_mult / total_stats.time_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - acceptance test time      : " << total_stats.time_acceptance << " (" << total_stats.time_acceptance / total_stats.time_estimate * 100 << "% of estimate time)\n\n";

        results.smart_total = total_stats.time_learn + total_stats.time_estimate;
        results.smart_total_ideal_search = total_stats.time_learn + total_stats.time_estimate - total_stats.time_search;
        results.smart_total_ideal_search_store = total_stats.time_learn + total_stats.time_estimate - total_stats.time_search - total_stats.time_store;

        std::cout << "total time (ideal search)               : " << total_stats.time_learn + total_stats.time_estimate - total_stats.time_search << "\n";
        std::cout << "total time (ideal store + ideal search) : " << total_stats.time_learn + total_stats.time_estimate - total_stats.time_search - total_stats.time_store << "\n\n";
        std::cout << "total time                              : " << results.smart_total << "\n\n";

        std::cout << std::setprecision(2)
              << 100.00 *
                 double(total_stats.learning_counter) /
                 double(total_stats.total_counter ) << "% of training : "
              << total_stats.learning_counter << " cells out of "
              << total_stats.total_counter << "\n\n";


    }
    else{
        results.conv_total = total_stats.time_learn + total_stats.time_estimate;

        std::cout << std::setprecision(6)
                  <<  "total time                     : "
                  << results.conv_total << "\n\n";
    }

    // The output stream of the data file
    std::ofstream datafile;

    // The floating-point precision in the output.
    int precision = 6;

    // Document times
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    datafile.open(folder + "/" + file + "-steps.txt", std::ofstream::out | std::ofstream::trunc);
    // Output statuses such that steps' go vertically and cells horizontally
    if (datafile.is_open()) {
        /*
        if (smart) {
            datafile
                    << "time_transport   time_equilib    time_estimate   time_search     "
                    << "time_mat_vect   time_acceptance time_learn      time_store      "
                    << "time_gibbs_min tree   smart_counter\n";
        }else{
            datafile
                    << "time_transport   time_equilib    time_learn\n";
        }
        */
        datafile << std::scientific << std::setprecision(precision);
        unsigned length = rt_times.size();
        for (unsigned i = 0; i < length; i++){
            datafile << rt_times[i] << "\t ";
            datafile << eq_times[i] << "\t ";
            if (smart) {
                // Output statuses collected while stepping with ReactiveTransportSolver
                datafile << step_stats[i].time_estimate << "\t "
                          << step_stats[i].time_search << "\t "
                          << step_stats[i].time_mat_vect_mult << "\t "
                          << step_stats[i].time_acceptance << "\t "
                          << step_stats[i].time_learn << "\t "
                          << step_stats[i].time_store << "\t "
                          << step_stats[i].time_gibbs_min << "\t"
                          << step_stats[i].tree_size << "\t"
                          << step_stats[i].learning_counter << "\n";

            } else
                datafile << step_stats[i].time_learn << "\n";
        }

    }

    // Close the data file
    datafile.close();

}

/// Constructor for ReactiveTransportSolver class
ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem &system)
: system_(system)
{
    // Set boundary condition
    setBoundaryState(ChemicalState(system));
}

/// Set options of the equilibrium solver
auto ReactiveTransportSolver::setEquilibriumOptions(const EquilibriumOptions &opts) -> void {

    if (options.smart)  {
        Assert(smart_equilibriumsolver != nullptr,
                "Smart Equilibrium Solver has not been initilized.",
                "Call initialize() method before setting the options of the equilibrium solver.");
        smart_equilibriumsolver->setOptions(opts);}
    else{
        Assert(equilibriumsolver != nullptr,
               "Equilibrium Solver has not been initilized.",
               "Call initialize() method before setting the options of the equilibrium solver.");
        equilibriumsolver->setOptions(opts);
    }
}

auto ReactiveTransportSolver::setMesh(const Mesh &mesh) -> void {
    transportsolver.setMesh(mesh);
}

auto ReactiveTransportSolver::setVelocity(double val) -> void {
    transportsolver.setVelocity(val);
}

auto ReactiveTransportSolver::setDiffusionCoeff(double val) -> void {
    transportsolver.setDiffusionCoeff(val);
}

auto ReactiveTransportSolver::setBoundaryState(
        const ChemicalState &state) -> void {
    bbc = state.elementAmounts();
}

auto ReactiveTransportSolver::setTimeStep(double val) -> void {
    transportsolver.setTimeStep(val);
}

auto ReactiveTransportSolver::output() -> ChemicalOutput {
    outputs.emplace_back(ChemicalOutput(system_));
    return outputs.back();
}

auto ReactiveTransportSolver::initialize() -> void {

    // Initialize mesh and corresponding amount e
    const Mesh &mesh = transportsolver.mesh();
    const Index num_elements = system_.numElements();
    const Index num_cells = mesh.numCells();

    // Initialize amount of elements in fluid and solid phases
    bf.resize(num_cells, num_elements);
    bs.resize(num_cells, num_elements);
    b.resize(num_cells, num_elements);

    // Initialize equilibrium solver based on the parameter
    if (options.smart)  smart_equilibriumsolver = std::make_unique<SmartEquilibriumSolver>(system_);
    else        equilibriumsolver = std::make_unique<EquilibriumSolver>(system_);

    // Initialize equilibrium solver based on the parameter
    transportsolver.setOptions();
    transportsolver.initialize();
}

auto ReactiveTransportSolver::step(ChemicalField &field, ReactiveTransportResult& rt_result) -> void {


    const auto &mesh = transportsolver.mesh();
    const auto &num_elements = system_.numElements();
    const auto &num_cells = mesh.numCells();
    const auto &ifs = system_.indicesFluidSpecies();
    const auto &iss = system_.indicesSolidSpecies();

    Time start;
    /*
    if (options.smart){
        std::cout << "learning time                  : " << rt_result.equilibrium.smart.learn_stats.time_learn << std:: endl;
        std::cout << "   - gibbs minimization time   : " << rt_result.equilibrium.smart.learn_stats.time_gibbs_min << std:: endl;
        std::cout << "   - store time                : " << rt_result.equilibrium.smart.learn_stats.time_store << std:: endl;
        std::cout << "estimating time                : " << rt_result.equilibrium.smart.estimate_stats.time_estimate << std:: endl;
        std::cout << "   - search time               : " << rt_result.equilibrium.smart.estimate_stats.time_search << std:: endl;
        std::cout << "   - matrix-vector mult. time  : " << rt_result.equilibrium.smart.estimate_stats.time_mat_vect_mult << std:: endl;
        std::cout << "   - acceptance test time      : " << rt_result.equilibrium.smart.estimate_stats.time_acceptance << std:: endl;
        std::cout << "total equilibrium time         : " << rt_result.equilibrium.smart.learn_stats.time_learn
                                                            + rt_result.equilibrium.smart.estimate_stats.time_estimate << std:: endl << std:: endl;
    }
    else{
        std::cout << "learning time                  : " << std::setprecision(6) << rt_result.equilibrium.stats.time_learn << std:: endl;
    }
    */
    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Left boundary condition cell
    unsigned int icell_bc = 0;
    double phi_bc = field[icell_bc].properties().fluidVolume().val;

    // Start profiling reactive transport
    if (options.profiling) start = time();

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {

        // Scale BC with a porosity of the boundary cell
        transportsolver.setBoundaryValue(phi_bc * bbc[ielement]);
        transportsolver.step(bf.col(ielement));
    }
    // Sum the amounts of elements distributed among fluid and solid species
    b.noalias() = bf + bs;

    // End profiling for the reactive transport
    if (options.profiling)  rt_result.rt_time = elapsed(start);

    // Open the the file for outputting chemical states
    for (auto output : outputs) {
        output.suffix("-" + std::to_string(steps));
        output.open();
    }

    for (Index icell = 0; icell < num_cells; ++icell) {

        const double T = field[icell].temperature();
        const double P = field[icell].pressure();

        // Start profiling equlibrium
        if (options.profiling) start = time();

        if (options.smart) {
            // Solve with a smart equilibrium solver
            rt_result.equilibrium += smart_equilibriumsolver->solve(field[icell], T, P, b.row(icell), steps, icell);

            // End profiling for the equilibrium calculations (accumulate cell-wise)
            if (options.profiling)  rt_result.eq_time += elapsed(start);

            // Update the time spend for either for learning or estimating
            if (!rt_result.equilibrium.smart.succeeded){
                rt_result.equilibrium.smart.addLearningIndex(icell);
                rt_result.equilibrium.smart.tree_size++;
            }
        } else {
            // Solve with a conventional equilibrium solver
            rt_result.equilibrium += equilibriumsolver->solve(field[icell], T, P, b.row(icell));

            // End profiling for the conventional equilibrium calculations (accumulate cell-wise)
            if (options.profiling)  rt_result.eq_time += elapsed(start);

        }

        for (auto output : outputs)
            output.update(field[icell], icell);
    }

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    if (options.smart && steps == rt_result.nsteps - 1)    smart_equilibriumsolver->showTree(steps);

    /*
    if (options.smart){
        for (auto indx : rt_result.equilibrium.smart.learning_states_indx)
            std::cout << indx << "\t";

        std::cout << std::endl << std::endl;

        std::cout << "learning time                  : " << rt_result.equilibrium.smart.learn_stats.time_learn << std:: endl;
        std::cout << "   - gibbs minimization time   : " << rt_result.equilibrium.smart.learn_stats.time_gibbs_min << std:: endl;
        std::cout << "   - store time                : " << rt_result.equilibrium.smart.learn_stats.time_store << std:: endl;
        std::cout << "estimating time                : " << rt_result.equilibrium.smart.estimate_stats.time_estimate << std:: endl;
        std::cout << "   - search time               : " << rt_result.equilibrium.smart.estimate_stats.time_search << std:: endl;
        std::cout << "   - matrix-vector mult. time  : " << rt_result.equilibrium.smart.estimate_stats.time_mat_vect_mult << std:: endl;
        std::cout << "   - acceptance test time      : " << rt_result.equilibrium.smart.estimate_stats.time_acceptance << std:: endl;
        std::cout << "total equilibrium time         : " << rt_result.equilibrium.smart.learn_stats.time_learn
                                                            + rt_result.equilibrium.smart.estimate_stats.time_estimate << std:: endl << std:: endl;
    }
    else{
        std::cout << "learning time                  : " << std::setprecision(6) << rt_result.equilibrium.stats.time_learn << std:: endl;
    }
    std::cout << "equilibrium time               : " << std::setprecision(6) << rt_result.eq_time << std:: endl;
    std::cout << "transport time                 : " << std::setprecision(6) << rt_result.rt_time << std:: endl << std:: endl;
    */
    // Collect the amounts of elements in the solid and fluid species
    ++steps;

}


} // namespace Reaktoro
