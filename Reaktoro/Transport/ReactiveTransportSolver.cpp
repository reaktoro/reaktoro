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
#include <numeric>
#include <tuple>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

/// Implementation of class EquilibriumProfiler
///
/// Constructor for the EquilibriumProfiler class
EquilibriumProfiler::EquilibriumProfiler(Profiling what): Profiler(what){}

/// Update vector with learning times
auto EquilibriumProfiler::updateLearning(int step) -> void
{
    /// Stop profiling
    std::chrono::duration<double> elapsed = clock::now() - start;
    double learn_time = elapsed.count();

    /// If the vector already has been initialized (on the current step), add the time to the existing value
    /// otherwise, emplace back the value to the vector as a new time on the current step
    if (!learn_times.empty() && learn_times.size() > step)
        learn_times.at(step) += learn_time;
    else
        learn_times.emplace_back(learn_time);

    /// If the vector with estimate_times has not been initialized at the current step
    /// add the zeros vector
    if (estimate_times.size() == step)
        estimate_times.emplace_back(std::vector<double>(estimate_data, 0.0));

    /// Increase the size of the current tree (storing reference states)
    tree_height ++;
}

/// Update vector with estimating times
auto EquilibriumProfiler::updateEstimating(int step, std::vector<double> time_partition) -> void{

    /// Stop profiling
    std::chrono::duration<double> elapsed =clock::now() - start;
    double est_time = elapsed.count();

    /// If the vector with estimate_times has not been initialized at the current step add the zeros vector
    if (estimate_times.size() == step) estimate_times.push_back(std::vector<double>(4, 0.0));
    /// Update statistics
    estimate_times.at(step)[0] += est_time;             // estimate time
    estimate_times.at(step)[1] += time_partition[0];    // search time
    estimate_times.at(step)[2] += time_partition[1];    // data fetch time
    estimate_times.at(step)[3] += time_partition[2];    // matrix-vector time
    estimate_times.at(step)[4] = tree_height;           // reference states tree

    /// If the vector with learn_times has not been initialized at the current step add zero to the vector
    if (learn_times.size() == step) learn_times.emplace_back(0.0);
}

/// Summary output to the console
auto EquilibriumProfiler::consoleOutput(int step) -> void
{
    // 0 - time for estimating
    // 1 - time for ref.element search
    // 2 - time for data fetching
    // 3 - time for matrix-vector manipulation
    std::cout << "total estimate time : " << estimate_times[step][0] << ""
              << " - ref.element search : " << estimate_times[step][1] / estimate_times[step][0] * 100 << "% "
              << " - data fetching : " << estimate_times[step][2] / estimate_times[step][0] * 100 << "% "
              << " - matrix-vector oper. : " << estimate_times[step][3] / estimate_times[step][0] * 100 << "% \n";
    std::cout << "total learning time : " << learn_times[step] << " \n\n";
}

/// Output the profiling results to the file
auto EquilibriumProfiler::fileOutput(const std::string & file) -> void
{
    /// The output stream of the data file
    std::ofstream datafile;

    /// The floating-point precision in the output.
    int precision = 6;

    /// Open the data file
    if(!file.empty())
        datafile.open(file + "-EQ-CW.txt", std::ofstream::out | std::ofstream::trunc);

    /// Output the header of the data file
    if(datafile.is_open())
    {
        /// Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);

        /// Output times collected while profiling
        unsigned int length = learn_times.size();
        for (unsigned int i = 0; i < length; i++)
            datafile << learn_times[i] + estimate_times[i][0] << "\t "
                      << learn_times[i] + estimate_times[i][0] - estimate_times[i][1] << "\n";
    }
    
    datafile.close();

    /// Open the date file for the estimation time analysis
    if(!file.empty())
        datafile.open(file + "-EQ-EST.txt", std::ofstream::out | std::ofstream::trunc);

    /// Output the header of the data file
    if(datafile.is_open())
    {
        /// Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);
        /// Output times collected while profiling
        unsigned int length = learn_times.size();
        for (unsigned int i = 0; i < length; i++)
            /// Output estimate time, search time, and the tree height
            datafile << estimate_times[i][0] << "\t "
                     << estimate_times[i][1] << "\t "
                     << estimate_times[i][4] << "\n";
    }

    datafile.close();

    /// Open the date file for the estimation time analysis
    if(!file.empty())
        datafile.open(file + "-EQ-MV.txt", std::ofstream::out | std::ofstream::trunc);

    /// Output the header of the data file
    if(datafile.is_open())
    {
        /// Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);

        /// Output times collected while profiling
        unsigned int length = learn_times.size();
        for (unsigned int i = 0; i < length; i++)
            /// Output estimate time, search time, and the tree height
            datafile << estimate_times[i][2] << "\t " << estimate_times[i][3] << "\n";
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
    /// The output filestream
    std::ofstream datafile;

    /// The suffix of the datafile
    std::string suffix;

    /// The floating-point precision in the output.
    int precision = 6;

    /// The suffix of the output file is dependent on Profiling
    switch (this->getProfilingSubject())
    {
        case Profiling::RT:     suffix = "RT"; break;
        case Profiling::EQ:     suffix = "EQ"; break;
        case Profiling::CK:     suffix = "CK";  break;
        case Profiling::Total:  suffix = "Total"; break;
        case Profiling::EQ_CW:  suffix = "EQ-CW"; break;
    }
    /// Open the data file
    if(!file.empty())
        datafile.open(file + "-" + suffix + ".txt",
                      std::ofstream::out | std::ofstream::trunc);
    /// Output the header of the data file
    if(datafile.is_open()) {
        /// Set scientific mode and defined precision
        datafile << std::scientific << std::setprecision(precision);
        /// Output times collected while profiling
        for (double time : times)   datafile << time << "\n";
    }
    /// Close output filestream
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
SolverStatus::SolverStatus(const std::string & folder, const std::string & file) 
: folder(folder), file(file)
{}

/// Output statuses of all steps and all cells to a file
auto SolverStatus::output(const Index & i) -> void
{
    /// The output stream of the data file
    std::ofstream datafile;

    /// Statuses for opening of the file
    auto opt = std::ofstream::out;

    /// Depending on the step, open new file or append to existing one
    if(i != 0 && !file.empty())         opt |= std::ofstream::app;
    else if(i == 0 && !file.empty())    opt |= std::ofstream::trunc;
    datafile.open(folder + "/" + file + ".txt", opt);

    /// Output statuses such that steps' go vertically and cells horizontally
    if(datafile.is_open()) {
        /// Output statuses collected while stepping with ReactiveTransportSolver
        for (bool est : statuses) datafile << std::to_string(est) << "\t";
        datafile << "\n";
    }
    
    /// Clear collected statuses
    statuses.clear();

    /// Close the data file
    datafile.close();
}
/// Constructor for ReactiveTransportSolver class
ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem &system, const bool &is_smart)
: system_(system), smart(is_smart)
{
    /// Define equilibrium solver based on the parameter
    if (smart)  smart_equilibriumsolver = SmartEquilibriumSolver(system);
    else        equilibriumsolver = EquilibriumSolver(system);

    /// Set boundary condition
    setBoundaryState(ChemicalState(system));
}

/// Set options of the Equlilibrium solver
auto ReactiveTransportSolver::setEquilibriumOptions(const EquilibriumOptions &options) -> void {
    if (smart)  smart_equilibriumsolver.setOptions(options);
    else        equilibriumsolver.setOptions(options);
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

auto ReactiveTransportSolver::cellprofile(Profiling subject) -> EquilibriumProfiler {
    eq_cell_profiler = std::make_unique<EquilibriumProfiler>(subject);
    return *eq_cell_profiler;
}
auto ReactiveTransportSolver::profile(Profiling subject) -> Profiler {
    profilers.push_back(Profiler(subject));
    return profilers.back();
}

auto ReactiveTransportSolver::trackStatus(const std::string &folder,
                                          const std::string &file) -> SolverStatus {
    status_trackers.push_back(SolverStatus(folder, file));
    return status_trackers.back();
}

auto ReactiveTransportSolver::outputProfiling(const std::string &folder) -> void {
    auto eq_profiler = find(begin(profilers), end(profilers), Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->fileOutput(folder);

    auto rt_profiler = std::find(begin(profilers), end(profilers), Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->fileOutput(folder);

    if (eq_cell_profiler) eq_cell_profiler->fileOutput(folder);
}

auto ReactiveTransportSolver::outputProfiling() -> void
{
    unsigned length = eq_cell_profiler->learn_times.size();

    double learn_time = 0;
    double matrix_vector_time = 0;
    double data_retreaval_time = 0;
    double minsearch_time = 0;
    double estimate_time = 0 ;

    for (unsigned i = 0; i < length; i++)
    {
        learn_time += eq_cell_profiler->learn_times[i];
        minsearch_time += eq_cell_profiler->estimate_times[i][1];
        estimate_time += eq_cell_profiler->estimate_times[i][0];
        data_retreaval_time += eq_cell_profiler->estimate_times[i][2];
        matrix_vector_time += eq_cell_profiler->estimate_times[i][3];
    }

    if (smart) {
        std::cout << std::setprecision(4);
        std::cout << "learning time                  : " << learn_time << " (" << learn_time / (learn_time + estimate_time) * 100 << "% of total time)\n";
        std::cout << "estimating time                : " << estimate_time  << " (" << estimate_time / (learn_time + estimate_time) * 100 << "% of total time)\n";
        std::cout << "   - search time               : " << minsearch_time << " (" << minsearch_time / estimate_time * 100 << "% of estimate time)\n";
        std::cout << "   - matrix-vector mult. time  : " << matrix_vector_time << " (" << matrix_vector_time / estimate_time * 100 << "% of estimate time)\n\n";

        std::cout << "total time                     : " << learn_time + estimate_time << "\n\n";

        std::cout << "estimating time (ideal search) : " << estimate_time - minsearch_time << "\n";
        std::cout << "total time (ideal search)      : " << learn_time + estimate_time - minsearch_time << "\n\n";
        std::cout << std::setprecision(2)
                  << 100.00 *
                     double(status_trackers[0].total_counter - status_trackers[0].smart_counter) /
                     double(status_trackers[0].total_counter) << "% of training : "
                  << status_trackers[0].total_counter - status_trackers[0].smart_counter << " cells out of "
                  << status_trackers[0].total_counter << "\n\n";

    }
    else
        std::cout << "total time                     : " << learn_time + estimate_time << "\n\n";

}

auto ReactiveTransportSolver::initialize() -> void {

    const Mesh &mesh = transportsolver.mesh();
    const Index num_elements = system_.numElements();
    const Index num_cells = mesh.numCells();

    bf.resize(num_cells, num_elements);
    bs.resize(num_cells, num_elements);
    b.resize(num_cells, num_elements);

    //transportsolver.initialize();
    transportsolver.initializeFullImplicit();
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

auto ReactiveTransportSolver::step_tracked(ChemicalField &field) -> void {
    const auto &mesh = transportsolver.mesh();
    const auto &num_elements = system_.numElements();
    const auto &num_cells = mesh.numCells();
    const auto &ifs = system_.indicesFluidSpecies();
    const auto &iss = system_.indicesSolidSpecies();

    //if (steps) std::cout << "system_ = " << system_ << std::endl;

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Porosity in the boundary cell

    std::vector<double> phi(num_cells, 0.0);
    for (Index icell = 0; icell < num_cells; ++icell) {
        phi[icell] = field[icell].properties().fluidVolume().val;
    }
    /*
    std::cout << "phi = \n";
    for (auto phi_i : phi)
        std::cout << phi_i << "\t";
    std::cout << std::endl;
    */

    unsigned int icell_bc = 0;
    double phi_bc = field[icell_bc].properties().fluidVolume().val;

    //std::cout << "field[icell].properties().volume().val = " << field[0].properties().volume().val << std::endl;

    // Find profiler for the reactive transort and start profiling
    auto rt_profiler = std::find(begin(profilers), end(profilers), Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->startProfiling();

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {

        // transportsolver.setBoundaryValue(bbc[ielement]);
        // Scale BC with a porousity of the boundary cell
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

        //std::cout << "b.row(icell) = \n" << b.row(icell) << std::endl;
        //std::cout << "field[icell] = \n" << field[icell] << std::endl;

        eq_cell_profiler->startProfiling();

        if (smart) {
            // Solve with a smart equilibrium solver
            EquilibriumResult res = smart_equilibriumsolver.solve(field[icell],
                                                                  T, P,
                                                                  b.row(icell));

            // Save the statuses of the cells depending on whether smart estimation or learning was triggered
            if (!status_trackers.empty())
                status_trackers[0].statuses.emplace_back(res.smart.succeeded);
            // Update the time spend for either for learning or estimating
            if (res.smart.succeeded) {
                eq_cell_profiler->updateEstimating(steps, res.smart.times);
                status_trackers[0].smart_counter++;
            }
            else
                eq_cell_profiler->updateLearning(steps);

            /*
            std::cout << "steps : " << steps << ", icell : " << icell  << " smart.succeeded : " << res.smart.succeeded;
            if (res.smart.succeeded)
                std::cout << ", estimate (s) : " << eq_cell_profiler->estimate_times[steps][0]
                          << ", tree height : " << eq_cell_profiler->tree_height << "\n";
            else
                std::cout << ", learning (s) : " << eq_cell_profiler->learn_times[steps] << "\n";
            */
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

    //eq_cell_profiler->consoleOutput(steps);

    //std::cout << "field[num_cells - 1] = \n" << field[num_cells - 1] << std::endl;

    // End profiling for the equllibrium calculations
    if (eq_profiler != end(profilers)) eq_profiler->endProfiling();

    // Output the states of the smart equilibrium algorithm to the file
    if (smart && !status_trackers.empty()) status_trackers[0].output(steps);

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    // Collect the amounts of elements in the solid and fluid species
    /*
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        std::cout << icell << " icell = " << std::scientific << bf.row(icell) << std::defaultfloat << std::endl;
    }
    */
    ++steps;
}

} // namespace Reaktoro
