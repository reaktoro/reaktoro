#include "ReactiveTransportSolver.hpp"

// Reaktoro's includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Common/Exception.hpp>

// C++ includes
#include <algorithm>
#include <fstream>
#include <iomanip>

namespace Reaktoro {

// Implementation of class Profiler
Profiler::Profiler(Reaktoro::Profiling subject_)
        : subject(subject_){
}
auto Profiler::startProfiling()  -> void
{
    start = std::chrono::high_resolution_clock::now();
}
auto Profiler::endProfiling() -> void
{
    finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

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

ReactiveTransportSolver::ReactiveTransportSolver(
        const ChemicalSystem &system, const bool &is_smart)
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

    // Collect the amounts of elements in the solid and fluid species
    for (Index icell = 0; icell < num_cells; ++icell) {
        bf.row(icell) = field[icell].elementAmountsInSpecies(ifs);
        bs.row(icell) = field[icell].elementAmountsInSpecies(iss);
    }

    // Find profiler for the reactive transort and start profiling
    auto rt_profiler = std::find(begin(profilers), end(profilers),
                                 Profiling::RT);
    if (rt_profiler != end(profilers)) rt_profiler->startProfiling();

    // Transport the elements in the fluid species
    for (Index ielement = 0; ielement < num_elements; ++ielement) {
        transportsolver.setBoundaryValue(bbc[ielement]);
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
    auto eq_profiler = find(begin(profilers), end(profilers),
                            Profiling::EQ);
    if (eq_profiler != end(profilers)) eq_profiler->startProfiling();


    for (Index icell = 0; icell < num_cells; ++icell) {
        const double T = field[icell].temperature();
        const double P = field[icell].pressure();
        /*
        using VectorConstRef   = Eigen::Ref<const Eigen::VectorXd>;           /// < Alias to Eigen type Ref<const VectorXd>.
        std::cout << "b(icell) = \n";
        VectorConstRef vector = b.row(icell);
        int len = vector.size();
        std::cout << std::scientific;
        for (auto i = 0; i < len; i++)
            std::cout << vector[i] << "\n";
        std::cout << std::endl;

        std::cout << " field[icell] = \n";
        std::cout << field[icell] << std::endl;
        */
        if (smart) {
            // Solve with a smart equilibrium solver
            EquilibriumResult res = smart_equilibriumsolver.solve(
                    field[icell], T, P, b.row(icell));
            // Save the states of the cells depending on whether smart estimation or learning was triggered
            if (!status_trackers.empty())
                status_trackers[0].statuses.emplace_back(
                        res.smart.succeeded);
        } else {
            // Solve with a conventional equilibrium solver
            equilibriumsolver.solve(field[icell], T, P, b.row(icell));
        }
        for (auto output : outputs)
            output.update(field[icell], icell);
    }
    // End profiling for the equllibrium calculations
    if (eq_profiler != end(profilers)) eq_profiler->endProfiling();

    // Output the states of the smart equilibrium algorithm to the file
    if (smart && !status_trackers.empty()) status_trackers[0].output(steps);

    // Output chemical states in the output files
    for (auto output : outputs)
        output.close();

    ++steps;
}
}
