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
ReactiveTransportResult::ReactiveTransportResult(const int & ncells, const int & nsteps, const bool& smart ) : ncells(ncells), nsteps(nsteps), is_smart(smart){}
auto ReactiveTransportResult::resetTime() -> void
{

    if (is_smart) {
        // Reset the CPU times
        equilibrium.smart.estimate_stats.time_estimate = 0.0;
        equilibrium.smart.estimate_stats.time_search = 0;
        equilibrium.smart.estimate_stats.time_mat_vect_mult = 0;
        equilibrium.smart.estimate_stats.time_acceptance = 0;
        equilibrium.smart.learn_stats.time_learn = 0;
        equilibrium.smart.learn_stats.time_store = 0;
        equilibrium.smart.learn_stats.time_gibbs_min = 0;
        equilibrium.smart.tree_size = 0;
        equilibrium.smart.smart_counter = 0;

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
ReactiveTransportProfiler::ReactiveTransportProfiler(const std::string &folder, const bool & smart)
    : folder(folder), smart(smart) {}

/// Process results collected per each step of reactive transort
/// \param rt_class that stores the results of the reactive transport
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
        stats.smart_counter = smart_res.smart_counter;
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
        total_stats.smart_counter += smart_res.smart_counter;
        total_stats.total_counter += rt_result.ncells;

        // Fill in the statuses array based on the
        auto it = smart_res.learning_states_indx.begin();

        for (Index i = 0; i < rt_result.ncells; i++) {
            if (it != smart_res.learning_states_indx.end() && i == *it)
                statuses.emplace_back(false);
            else
                statuses.emplace_back(true);
        }

        // Null the CUP time
        rt_result.resetTime();

    }
    else{
        EquilibriumResult res = rt_result.equilibrium;
        Statistics stats;

        // Initialize the results step-wise
        stats.time_learn = res.stats.time_learn;
        step_stats.emplace_back(stats);

        rt_result.resetTime();
        // Accumulate the step-wise results
        total_stats.time_learn += res.stats.time_learn;
    }
    // Add reactive transport and equilibrium time
    rt_times.emplace_back(rt_result.rt_time);
    eq_times.emplace_back(rt_result.eq_time);
};

/// Output results collected per each step of reactive transport
/// \param file_name of the file to be updated
auto ReactiveTransportProfiler::output(const std::string &file,
                                       const Index &step) -> void {
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
        datafile.open(folder + "/" + file + "-statuses.txt", opt);
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

    // Document times
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    datafile.open(folder + "/" + file + "-steps.txt", opt);
    // Output statuses such that steps' go vertically and cells horizontally
    if (datafile.is_open()) {
        if (smart) {
            // Output statuses collected while stepping with ReactiveTransportSolver
            Statistics stats = step_stats.back();
            datafile << std::scientific << std::setprecision(precision);
            datafile << stats.time_estimate << "\t "
                     << stats.time_search << "\t "
                     << stats.time_mat_vect_mult << "\t "
                     << stats.time_acceptance << "\t "
                     << stats.time_learn << "\t "
                     << stats.time_store << "\t "
                     << stats.time_gibbs_min << "\t"
                     << stats.tree_size << "\t"
                     << stats.smart_counter << "\n";

        } else{
            Statistics stats = step_stats.back();
            datafile << stats.time_learn << "\n";
        }
    }

    // Close the data file
    datafile.close();

};

auto ReactiveTransportProfiler::console(const Index & step) -> void
{
    if (smart){
        std::cout << "total estimate time : " << step_stats[step].time_estimate << ": "
        << " - ref.element search : " << step_stats[step].time_search / step_stats[step].time_estimate * 100 << "%, "
        << " - matrix-vector oper. : " << step_stats[step].time_mat_vect_mult / step_stats[step].time_estimate * 100 << "%, "
        << " - acceptance. : " << step_stats[step].time_acceptance / step_stats[step].time_estimate * 100 << "% \n";
        std::cout << "total learning time : " << step_stats[step].time_learn << ": "
        << " - gibbs min. : " << step_stats[step].time_gibbs_min / step_stats[step].time_learn * 100 << "%, "
        << " - ref.element store : " << step_stats[step].time_store / step_stats[step].time_learn * 100 << "%\n";
    } else{
        std::cout << "total learning time : " << step_stats[step].time_learn;
    }
}
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
                 double(total_stats.total_counter - total_stats.smart_counter) /
                 double(total_stats.total_counter ) << "% of training : "
              << total_stats.total_counter - total_stats.smart_counter << " cells out of "
              << total_stats.total_counter << "\n\n";


    }
    else{
        results.conv_total = total_stats.time_learn + total_stats.time_estimate;

        std::cout << std::setprecision(6)
                  <<  "total time                     : "
                  << results.conv_total << "\n\n";
    }

}
/*
/// Implementation of ...
///
auto ReactiveTransportResult::output(const std::string & file) -> void
{

};
*/

/// Implementation of class EquilibriumProfiler
///
/// Constructor for the EquilibriumProfiler class
//EquilibriumProfiler::EquilibriumProfiler(Profiling what): Profiler(what){}
/*
/// Update vector with learning times (when conventional method is used)
auto EquilibriumProfiler::updateLearning(int step) -> void {

    // Stop profiling
    std::chrono::duration<double> elapsed = clock::now() - start;
    double learn_time = elapsed.count();

    // If the vector was not initialized, emplace back statistics' zero values on the current step
    if (learning_stats.size() == step || learning_stats.empty()) learning_stats.emplace_back(Statistics());
    learning_stats.at(step).time_learn += learn_time;

}

/// Update vector with learning times
auto EquilibriumProfiler::updateLearning(int step, SmartEquilibriumResult res) -> void {

    /// Stop profiling
    std::chrono::duration<double> elapsed = clock::now() - start;
    double learn_time = elapsed.count();

    /// If the vector was not initialized, emplace back statistics' zero values on the current step
    if (learning_stats.size() == step || learning_stats.empty()) learning_stats.emplace_back(Statistics());

    /// If the vector with estimate_times has not been initialized at the current step
    /// add the zeros vector
    if (estimating_stats.size() == step) estimating_stats.emplace_back(Statistics());

    /// Increase the size of the current tree (storing reference states)
    tree_height ++;
}


/// Update vector with estimating times
auto EquilibriumProfiler::updateEstimating(int step, SmartEquilibriumResult::EstimateStatistics stats) -> void{

    /// Stop profiling
    std::chrono::duration<double> elapsed =clock::now() - start;
    double est_time = elapsed.count();

    /// If the vector with estimate_times has not been initialized at the current step add the zeros vector
    if (estimating_stats.size() == step || estimating_stats.empty()) estimating_stats.emplace_back(Statistics());

    /// Update statistics
    estimating_stats.at(step).time_estimate += est_time;
    estimating_stats.at(step).time_search += stats.time_search;
    estimating_stats.at(step).time_matrix_vector_mult += stats.time_matrix_vector_mult;
    estimating_stats.at(step).time_acceptance_test += stats.time_acceptance_test;
    estimating_stats.at(step).tree_size = tree_height;

    /// If the vector with learn_times has not been initialized at the current step add zero to the vector
    if (learning_stats.size() == step) learning_stats.emplace_back(Statistics());
}
*/
/*
/// Summary output to the console
auto EquilibriumProfiler::consoleOutput(int step) -> void
{
    std::cout << "total estimate time : " << estimating_stats[step].time_estimate << ""
              << " - ref.element search : " << estimating_stats[step].time_search / estimating_stats[step].time_estimate * 100 << "% "
              << " - matrix-vector oper. : " << estimating_stats[step].time_matrix_vector_mult / estimating_stats[step].time_estimate * 100 << "% \n";
    std::cout << "total learning time : " << learning_stats[step].time_learn << " \n\n";
    std::cout << " - gibbs min. : " << learning_stats[step].time_gibbs_minimization << " \n\n";
    std::cout << " - ref.element store : " << learning_stats[step].time_store << " \n\n";
}

/// Output the profiling results to the file
auto EquilibriumProfiler::fileOutput(const std::string & file) -> void
{
    // The output stream of the data file
    std::ofstream datafile;

    // The floating-point precision in the output.
    int precision = 6;

    // Open the data file
    if(!file.empty())
        datafile.open(file + "-EQ-CW.txt", std::ofstream::out | std::ofstream::trunc);

    // Output the header of the data file
    if(datafile.is_open())
    {
        // Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);

        // Output times collected while profiling
        unsigned int length = learning_stats.size();
        for (unsigned int i = 0; i < length; i++)
            datafile << learning_stats[i].time_learn + estimating_stats[i].time_estimate << "\t "
                      << learning_stats[i].time_learn + estimating_stats[i].time_estimate - estimating_stats[i].time_search << "\n";
    }
    
    datafile.close();

    // Open the date file for the estimation time analysis
    if(!file.empty())
        datafile.open(file + "-EQ-EST.txt", std::ofstream::out | std::ofstream::trunc);

    // Output the header of the data file
    if(datafile.is_open())
    {
        // Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);
        // Output times collected while profiling
        unsigned int length = learning_stats.size();
        for (unsigned int i = 0; i < length; i++)
            // Output estimate time, search time, and the tree height
            datafile << estimating_stats[i].time_estimate << "\t "
                     << estimating_stats[i].time_search << "\t "
                     << estimating_stats[i].tree_size << "\n";
    }

    datafile.close();

}

/// Implementation of class EquilibriumProfiler
///
/// Constructors and destructors for the Profiler class
Profiler::Profiler(Reaktoro::Profiling subject_) : subject(subject_){}

/// Method to start profiling
auto Profiler::startProfiling() -> void
{
    start = clock::now();
}

/// Method to stop profiling
auto Profiler::endProfiling() -> void
{
    std::chrono::duration<double> elapsed = clock::now() - start;
    times.emplace_back(elapsed.count());
}

/// Method to output profiling results to the file
auto Profiler::fileOutput(const std::string & file) -> void
{
    // The output filestream
    std::ofstream datafile;

    // The suffix of the datafile
    std::string suffix;

    // The floating-point precision in the output.
    int precision = 6;

    // The suffix of the output file is dependent on Profiling
    switch (this->getProfilingSubject())
    {
        case Profiling::RT:     suffix = "RT"; break;
        case Profiling::EQ:     suffix = "EQ"; break;
        case Profiling::CK:     suffix = "CK";  break;
        case Profiling::Total:  suffix = "Total"; break;
        case Profiling::EQ_CW:  suffix = "EQ-CW"; break;
    }
    // Open the data file
    if(!file.empty())
        datafile.open(file + "-" + suffix + ".txt",
                      std::ofstream::out | std::ofstream::trunc);
    // Output the header of the data file
    if(datafile.is_open()) {
        // Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);
        // Output times collected while profiling
        for (double time : times)   datafile << time << "\n";
    }
    // Close output filestream
    datafile.close();
}

/// Get profiling subject
auto Profiler::getProfilingSubject() const -> Profiling
{
    return subject;
}

/// Operator == overwritten for fetching needed profiler from the vector profilers
auto Profiler::operator==(const Profiler& p) const -> bool
{
    return p.getProfilingSubject() == this->subject;
}


/// Implementation of class SolverStatus
///
/// Class for tracking the statuses of SmartEquilibriumSolver
SmartSolverStatus::SmartSolverStatus(const std::string & folder, const std::string & file)
: folder(folder), file(file)
{}

/// Output statuses of all steps and all cells to a file
auto SmartSolverStatus::output(const Index & i) -> void
{
    // The output stream of the data file
    std::ofstream datafile;

    // Statuses for opening of the file
    auto opt = std::ofstream::out;

    // Depending on the step, open new file or append to existing one
    if(i != 0 && !file.empty())         opt |= std::ofstream::app;
    else if(i == 0 && !file.empty())    opt |= std::ofstream::trunc;
    datafile.open(folder + "/" + file + ".txt", opt);

    // Output statuses such that steps' go vertically and cells horizontally
    if(datafile.is_open()) {
        // Output statuses collected while stepping with ReactiveTransportSolver
        for (bool est : statuses) datafile << std::to_string(est) << "\t";
        datafile << "\n";
    }
    
    // Clear collected statuses
    statuses.clear();

    // Close the data file
    datafile.close();
}

*/

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
/*
auto ReactiveTransportSolver::cellprofile(Profiling subject) -> EquilibriumProfiler {
    eq_cell_profiler = std::make_unique<EquilibriumProfiler>(subject);
    return *eq_cell_profiler;
}
auto ReactiveTransportSolver::profile(Profiling subject) -> Profiler {
    profilers.push_back(Profiler(subject));
    return profilers.back();
}

auto ReactiveTransportSolver::trackStatus(const std::string &folder,
                                          const std::string &file) -> SmartSolverStatus {
    status_trackers.push_back(SmartSolverStatus(folder, file));
    return status_trackers.back();
}

auto ReactiveTransportSolver::outputProfiling(const std::string &folder) -> void {
    auto eq_profiler = find(begin(profilers), end(profilers), Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->fileOutput(folder);

    auto rt_profiler = std::find(begin(profilers), end(profilers), Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->fileOutput(folder);

    if (eq_cell_profiler && options.smart) eq_cell_profiler->fileOutput(folder);
}

auto ReactiveTransportSolver::outputProfiling() -> void
{
    unsigned length = eq_cell_profiler->learning_stats.size();

    double total_estimate(0.0), total_learn(0.0);
    double total_matrix_vector_mult(0.0), total_search(0.0), total_acceptance_test(0.0);
    double total_store(0.0), total_gibbs_min(0.0);

    for (unsigned i = 0; i < length; i++) {
        total_learn += eq_cell_profiler->learning_stats[i].time_learn;

        if (options.smart) {
            total_estimate += eq_cell_profiler->estimating_stats[i].time_estimate;
            total_search += eq_cell_profiler->estimating_stats[i].time_search;
            total_matrix_vector_mult += eq_cell_profiler->estimating_stats[i].time_matrix_vector_mult;
            total_acceptance_test += eq_cell_profiler->estimating_stats[i].time_acceptance_test;
            total_gibbs_min += eq_cell_profiler->learning_stats[i].time_gibbs_minimization;
            total_store += eq_cell_profiler->learning_stats[i].time_store;
        }
    }

    if (options.smart) {
        std::cout << std::setprecision(6);

        std::cout << "learning time                  : " << total_learn << " (" << total_learn / (total_learn + total_estimate) * 100 << "% of total time)\n";
        std::cout << "   - gibbs minimization time   : " << total_gibbs_min << " (" << total_gibbs_min / total_learn * 100 << "% of learning time)\n";
        std::cout << "   - store time                : " << total_store << " (" << total_store / total_learn * 100 << "% of learning time)\n";
        std::cout << "estimating time                : " << total_estimate  << " (" << total_estimate / (total_learn + total_estimate) * 100 << "% of total time)\n";
        std::cout << "   - search time               : " << total_search << " (" << total_search / total_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - matrix-vector mult. time  : " << total_matrix_vector_mult << " (" << total_matrix_vector_mult / total_estimate * 100 << "% of estimate time)\n";
        std::cout << "   - acceptance test time      : " << total_acceptance_test << " (" << total_acceptance_test / total_estimate * 100 << "% of estimate time)\n\n";

        std::cout << "total time                     : " << total_learn + total_estimate << "\n\n";

        std::cout << "estimating time (ideal search) : " << total_estimate - total_search << "\n";
        std::cout << "total time (ideal search)      : " << total_learn + total_estimate - total_search << "\n\n";

        std::cout << std::setprecision(2)
                  << 100.00 *
                     double(status_trackers[0].total_counter - status_trackers[0].smart_counter) /
                     double(status_trackers[0].total_counter) << "% of training : "
                  << status_trackers[0].total_counter - status_trackers[0].smart_counter << " cells out of "
                  << status_trackers[0].total_counter << "\n\n";

    }
    else
        std::cout << std::setprecision(6)
                  <<  "total time                     : " << total_learn + total_estimate << "\n\n";

}
*/
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

/*
auto ReactiveTransportSolver::step(ChemicalField &field) -> void {
    const auto &mesh = transportsolver.mesh();
    const auto &num_elements = system_.numElements();
    const auto &num_cells = mesh.numCells();
    const auto &ifs = system_.indicesFluidSpecies();
    const auto &iss = system_.indicesSolidSpecies();

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Porosity in the boundary cell
    unsigned int icell_bc = 0;
    double phi_bc = field[icell_bc].properties().fluidVolume().val;

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {
        transportsolver.setBoundaryValue(phi_bc * bbc[ielement]);
        transportsolver.step(bf.col(ielement));
    }
    // Sum the amounts of elements distributed among fluid and solid species
    b.noalias() = bf + bs;

    // Open the the file for outputting chemical states
    for (auto output : outputs) {
        output.suffix("-" + std::to_string(steps));
        output.open();
    }

    for (Index icell = 0; icell < num_cells; ++icell) {
        const double T = field[icell].temperature();
        const double P = field[icell].pressure();

        // Solve with a smart or conventional equilibrium solver
        equilibriumsolver.solve(field[icell], T, P, b.row(icell));

        for (auto output : outputs)
            output.update(field[icell], icell);
    }

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    ++steps;
}

auto ReactiveTransportSolver::step(ChemicalField &field) -> void {

    const auto &mesh = transportsolver.mesh();
    const auto &num_elements = system_.numElements();
    const auto &num_cells = mesh.numCells();
    const auto &ifs = system_.indicesFluidSpecies();
    const auto &iss = system_.indicesSolidSpecies();

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Left boundary condition cell
    unsigned int icell_bc = 0;
    double phi_bc = field[icell_bc].properties().fluidVolume().val;

    // Find profiler for the reactive transort and start profiling
    auto rt_profiler = std::find(begin(profilers), end(profilers), Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->startProfiling();

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {

        // Scale BC with a porosity of the boundary cell
        transportsolver.setBoundaryValue(phi_bc * bbc[ielement]);
        transportsolver.step(bf.col(ielement));
    }
    // Sum the amounts of elements distributed among fluid and solid species
    b.noalias() = bf + bs;

    // End profiling for the reactive transpor
    if (rt_profiler != end(profilers)) rt_profiler->endProfiling();

    // Open the the file for outputting chemical states
    for (auto output : outputs) {
        output.suffix("-" + std::to_string(steps));
        output.open();
    }

    // Find profiler for the chemical equilibrium and start profiling
    auto eq_profiler = find(begin(profilers), end(profilers),Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->startProfiling();

    for (Index icell = 0; icell < num_cells; ++icell) {

        const double T = field[icell].temperature();
        const double P = field[icell].pressure();

        eq_cell_profiler->startProfiling();

        if (smart) {
            // Solve with a smart equilibrium solver
            EquilibriumResult res = smart_equilibriumsolver.solve(field[icell], T, P, b.row(icell));

            // Save the statuses of the cells depending on whether smart estimation or learning was triggered
            if (!status_trackers.empty())
                status_trackers[0].statuses.emplace_back(res.smart.succeeded);
            // Update the time spend for either for learning or estimating
            if (res.smart.succeeded) {
                eq_cell_profiler->updateEstimating(steps, res.smart.estimate_stats);
                status_trackers[0].smart_counter++;
            }
            else
                eq_cell_profiler->updateLearning(steps, res.smart.learn_stats);

                //std::cout << "b.row(icell)        :\n" << b.row(icell) << std::endl;
                //std::cout << "field[icell]        :\n" << field[icell] << std::endl;
        } else {
            // Solve with a conventional equilibrium solver
            equilibriumsolver.solve(field[icell], T, P, b.row(icell));
            eq_cell_profiler->updateLearning(steps);
        }

        for (auto output : outputs)
            output.update(field[icell], icell);
    }
    
    // Counter all the statuses
    status_trackers[0].total_counter += num_cells;

    // End profiling for the equllibrium calculations
    if (eq_profiler != end(profilers)) eq_profiler->endProfiling();

    // Output the states of the smart equilibrium algorithm to the file
    if (smart && !status_trackers.empty()) status_trackers[0].output(steps);

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    // Collect the amounts of elements in the solid and fluid species
    ++steps;
}
*/
auto ReactiveTransportSolver::step(ChemicalField &field, ReactiveTransportResult& rt_result) -> void {


    const auto &mesh = transportsolver.mesh();
    const auto &num_elements = system_.numElements();
    const auto &num_cells = mesh.numCells();
    const auto &ifs = system_.indicesFluidSpecies();
    const auto &iss = system_.indicesSolidSpecies();

    Time start;

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
            if (rt_result.equilibrium.smart.succeeded)
                rt_result.equilibrium.smart.smart_counter ++;
            else{
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

    // Collect the amounts of elements in the solid and fluid species
    ++steps;

}


} // namespace Reaktoro
