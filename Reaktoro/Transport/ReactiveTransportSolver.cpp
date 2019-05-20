#include "ReactiveTransportSolver.hpp"

// Reaktoro's includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Common/Exception.hpp>

// C++ includes
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <tuple>

namespace Reaktoro {



// Implementation of class Profiler

EquilibriumProfiler::EquilibriumProfiler(Profiling what): Profiler(what){}

auto EquilibriumProfiler::updateLearning(int step) -> void
{
    std::chrono::duration<double> elapsed = clock::now() - start;
    if (!learn_times.empty() && learn_times.size() > step)
        learn_times.at(step) += elapsed.count();
    else
        learn_times.emplace_back(elapsed.count());

    tree_height ++;
}

auto EquilibriumProfiler::updateEstimating(int step, std::vector<double> time_partition) -> void{

    std::chrono::duration<double> elapsed = clock::now() - start;
    std::cout << elapsed.count() << "\t"
              << time_partition[0] << "\t"
              << time_partition[1] << "\t"
              << time_partition[2] << "\n";

    estimate_times.at(step)[0] += elapsed.count();
    estimate_times.at(step)[1] += time_partition[0];
    estimate_times.at(step)[2] += time_partition[1];
    estimate_times.at(step)[3] += time_partition[2];
}
auto EquilibriumProfiler::consoleOutput() -> void
{
    std::cout << "Step \t Learning \t Estimating";
    unsigned length = learn_times.size();
    for (unsigned i = 0; i < length; i++)
        std::cout << i + 1 << "\t" << learn_times[i] << "\t" << estimate_times[i][0];
}
auto EquilibriumProfiler::outputEstimateTime(int step) -> void
{
    if (!estimate_times.empty())
        // 1 - time for ref.element search
        // 2 - time for data fetching
        // 3 - time for matrix-vector manipulation
        std::cout << "out of total time " << estimate_times[step][0] << ": \n"
                  << "\t - ref.element search  : " << estimate_times[step][1] << "\n"
                  << "\t - data fetching       : " << estimate_times[step][2] << "\n"
                  << "\t - matrix-vector oper. : " << estimate_times[step][3] << "\n";
    else
        estimate_times.emplace_back(std::vector<double>(4, 0.0));
}

Profiler::Profiler(Reaktoro::Profiling subject_) : subject(subject_){}
auto Profiler::startProfiling()  -> void { start = clock::now(); }
auto Profiler::endProfiling() -> void
{
    std::chrono::duration<double> elapsed = clock::now() - start;

    times.emplace_back(elapsed.count());
}
auto Profiler::fileOutput(const std::string & file) -> void
{
    /// The output stream of the data file
    std::ofstream datafile;

    /// The suffix of the datafile
    std::string suffix;

    /// The flag that indicates if scientific format should be used.
    bool scientific = true;

    /// The floating-point precision in the output.
    int precision = 6;

    switch (this->getProfilingSubject())
    {
        case Profiling::RT:     suffix = "RT"; break;
        case Profiling::EQ:     suffix = "EQ"; break;
        case Profiling::CK:     suffix = "CK";  break;
        case Profiling::Total:  suffix = "Total"; break;
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
        for (double time : times)
            datafile << time << "\n";
    }

}
auto Profiler::consoleOutput() -> void
{
    std::for_each(begin(times), end(times),
                  [&](const double & value){ std::cout << value << "\t"; });
}
auto Profiler::getProfilingSubject() const -> Profiling {
    return subject;
}
// Operator == overwritten for fetching needed profiler from the vector profilers
auto Profiler::operator==(const Profiler& p) const -> bool{
    return p.getProfilingSubject() == this->subject;
}

// Class for tracking the statuses of SmartEquilibriumSolver
SolverStatus::SolverStatus(const std::string & folder, const std::string & file) :
        folder(folder), file(file) {}
// Output statuses of all steps and all cells to a file
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

    // Output statuses such that
    // vertically we have steps' number and horizontally we cells
    if(datafile.is_open()) {
        // Output statuses collected while stepping with RT
        for (bool est : statuses)   datafile << std::to_string(est) << "\t";
        datafile << "\n";
    }
    // Clear collected statuses
    statuses.clear();

    // Close the data file
    datafile.close();
}

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem &system, const bool &is_smart)
        : system_(system), smart(is_smart) {
    if (smart)  smart_equilibriumsolver = SmartEquilibriumSolver(system);
    else        equilibriumsolver = EquilibriumSolver(system);

    setBoundaryState(ChemicalState(system));
}

auto ReactiveTransportSolver::setEquilibriumOptions(
        const EquilibriumOptions &options) -> void {
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

auto ReactiveTransportSolver::outputProfiling(
        const std::string &folder) -> void {
    auto eq_profiler = find(begin(profilers), end(profilers), Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->fileOutput(folder);

    auto rt_profiler = std::find(begin(profilers), end(profilers), Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->fileOutput(folder);
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

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {

        transportsolver.setBoundaryValue(bbc[ielement]);
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

    if (steps) std::cout << "system_ = " << system_ << std::endl;

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        std::cout << "bf.row(icell) = " << bf.row(icell) << std::endl;
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Find profiler for the reactive transort and start profiling
    auto rt_profiler = std::find(begin(profilers), end(profilers),
                                 Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->startProfiling();

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {

        if (ielement == num_elements - 1) {
            std::cout << "bbc[ielement] = " << bbc[ielement] << std::endl;
            std::cout << "bf.col(ielement) = " << bf.col(ielement) << std::endl;
        }
        transportsolver.setBoundaryValue(bbc[ielement]);
        transportsolver.step(bf.col(ielement));

        if (ielement == num_elements - 1)
            std::cout << "bbf.col(ielement) = " << bf.col(ielement) << std::endl;

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
    auto eq_profiler = find(begin(profilers), end(profilers),
                            Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->startProfiling();

    for (Index icell = 0; icell < num_cells; ++icell) {
        const double T = field[icell].temperature();
        const double P = field[icell].pressure();

        std::cout << "b.row(icell) = \n" << b.row(icell) << std::endl;
        std::cout << "field[icell] = \n" << field[icell] << std::endl;

        if (smart) {
            eq_cell_profiler->startProfiling();
            // Solve with a smart equilibrium solver
            EquilibriumResult res = smart_equilibriumsolver.solve(field[icell],
                                                                  T, P,
                                                                  b.row(icell));

            // Save the statuses of the cells depending on whether smart estimation or learning was triggered
            if (!status_trackers.empty())
                status_trackers[0].statuses.emplace_back(res.smart.succeeded);

            // Update the time spend for either for learnign or estimating
            if (res.smart.succeeded)
                eq_cell_profiler->updateEstimating(steps, res.smart.times);
            else
                eq_cell_profiler->updateLearning(steps);

            std::cout << "icell               : " << icell << "\n";
            if (res.smart.succeeded) {
                std::cout << "ref.element search  : " << res.smart.times[0] << "\n"
                          << "data fetching       : " << res.smart.times[1] << "\n"
                          << "matrix-vector oper. : " << res.smart.times[2] << "\n";
            }
            else {
                std::cout << "learning time       : " << eq_cell_profiler->learn_times[steps] << "\n";
            }
            std::cout << "b.row(icell)        :\n" << b.row(icell) << std::endl;
            std::cout << "field[icell]        :\n" << field[icell] << std::endl;

        } else {
            // Solve with a conventional equilibrium solver
            equilibriumsolver.solve(field[icell], T, P, b.row(icell));
        }

        for (auto output : outputs)
            output.update(field[icell], icell);
    }
    eq_cell_profiler->outputEstimateTime(steps);

    //std::cout << "field[num_cells - 1] = \n" << field[num_cells - 1] << std::endl;

    // End profiling for the equllibrium calculations
    if (eq_profiler != end(profilers)) eq_profiler->endProfiling();

    // Output the states of the smart equilibrium algorithm to the file
    if (smart && !status_trackers.empty()) status_trackers[0].output(steps);

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        std::cout << "Z(icell) = " << bf.row(icell)[8] << std::endl;
    }
    ++steps;
}
}
