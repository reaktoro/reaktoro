// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "ChemicalPlot.hpp"

// C++ includes
#include <fstream>
#include <iomanip>

// Boost includes
#include <boost/format.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalQuantity.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

struct ChemicalPlot::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The reaction system instance
    ReactionSystem reactions;

    /// The chemical quantity instance
    ChemicalQuantity quantity;

    /// The name of the plot.
    std::string name;

    /// The quantity that spans the x-axis
    std::string x = "t";

    /// The quantities to be plotted along the y-axis
    std::vector<std::string> y = {"t"};

    /// The names of each curve given by member `y`.
    std::vector<std::string> legend;

    /// The Gnuplot commands used to configure the plot.
    std::string config;

    /// The frequency in which the plot is refreshed per second.
    unsigned frequency = 30;

    /// The name of the data file.
    std::string dataname;

    /// The name of the gnuplot script file.
    std::string plotname;

    /// The output stream of the data file.
    std::ofstream datafile;

    /// The output stream of the gnuplot script file.
    std::ofstream plotfile;

    /// The pointer to the pipe connecting to Gnuplot
    FILE* pipe = nullptr;

    Impl()
    {}

    Impl(const ChemicalSystem& system)
    : system(system), quantity(system)
    {}

    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions), quantity(reactions)
    {}

    ~Impl()
    {
        close();
    }

    auto open() -> void
    {
        // Ensure the plot is closed
        close();

        // Make sure name is not empty
        if(name.empty())
            name = "plot" + std::to_string(std::rand());

        // Make sure legend is not empty
        if(legend.empty())
            legend = y;

        // Initialize the names of the data and gnuplot script files
        dataname = name + ".dat";
        plotname = name + ".plt";

        // Open the data and gnuplot script files
        datafile.open(dataname);
        plotfile.open(plotname);

        // Output the name of each quantity in the data file
        datafile << std::left << std::setw(20) << x;
        for(auto yi : y)
            datafile << std::left << std::setw(20) << yi;
        datafile << std::endl;

        // Initialize the Gnuplot script file with the provided configuration
        plotfile << config;

        // Define a Gnuplot variable that is a list of titles for the curves
        plotfile << "titles = '";
        for(const auto& title : legend)
            plotfile << title << (&title != &legend.back() ? " " : "");
        plotfile << "'\n" << std::endl;

        // Define the formatted string that represents the plot part of the Gnuplot script
        std::string script =
            "previous = current\n"
            "current = system('%1% %2%')\n"
            "pause %3%\n"
            "plot for [i=2:%4%] '%2%' using 1:i with lines lt i-1 lw 2 title word(titles, i-1)\n"
            "if(current ne previous) reread";

        // On Windows, use the `dir` command on the data file to check its state.
        // On any other OS, use the `ls -l` command instead.
#if _WIN32
        std::string cmd = "dir";
#else
        std::string cmd = "ls -l";
#endif
        // Define auxiliary variables for the plot
        auto imax = 1 + y.size();
        auto wait = 1.0/frequency;

        // Finalize the Gnuplot script
        plotfile << boost::format(script) % cmd % dataname % wait % imax;

        // Flush the plot file to ensure its correct state before the plot starts
        plotfile.flush();
    }

    auto close() -> void
    {
        if(pipe != nullptr)
            pclose(pipe);
        pipe = nullptr;
    }

    auto update(const ChemicalState& state, double t) -> void
    {
        // Output the current chemical state to the data file.
        quantity.update(state, t);
        datafile << std::left << std::setw(20) << quantity.value(x);
        for(auto word : y)
            datafile << std::left << std::setw(20) << quantity.value(word);
        datafile << std::endl;

        // Open the Gnuplot plot after the first data has been output to the data file.
        // This ensures that Gnuplot opens the plot without errors/warnings.
        if(pipe == nullptr)
        {
            auto command = ("gnuplot -persist -e \"current=''\" " + plotname).c_str();
            pipe = popen(command, "w");
        }
    }
};

ChemicalPlot::ChemicalPlot()
: pimpl(new Impl())
{}

ChemicalPlot::ChemicalPlot(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalPlot::ChemicalPlot(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalPlot::~ChemicalPlot()
{}

auto ChemicalPlot::name(std::string name) -> void
{
    pimpl->name = name;
}

auto ChemicalPlot::x(std::string x) -> void
{
    pimpl->x = x;
}

auto ChemicalPlot::y(std::vector<std::string> y) -> void
{
    pimpl->y = y;
}

auto ChemicalPlot::y(std::string y) -> void
{
    pimpl->y = split(y);
}

auto ChemicalPlot::legend(std::vector<std::string> legend) -> void
{
    pimpl->legend = legend;
}

auto ChemicalPlot::legend(std::string legend) -> void
{
    pimpl->legend = split(legend);
}

auto ChemicalPlot::frequency(unsigned frequency) -> void
{
    pimpl->frequency = frequency;
}

auto ChemicalPlot::operator<<(std::string command) -> ChemicalPlot&
{
    pimpl->config.append(command + "\n");
    return *this;
}

auto ChemicalPlot::operator<<(std::stringstream command) -> ChemicalPlot&
{
    pimpl->config.append(command.str());
    return *this;
}

auto ChemicalPlot::open() -> void
{
    pimpl->open();
}

auto ChemicalPlot::update(const ChemicalState& state, double t) -> void
{
    pimpl->update(state, t);
}

} // namespace Reaktoro
