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
    ChemicalSystem system;

    ReactionSystem reactions;

    ChemicalQuantity quantity;

    ChemicalPlotOptions options;

    std::string dataname;

    std::string plotname;

    std::ofstream datafile;

    std::ofstream plotfile;

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

    auto open(const ChemicalPlotOptions& options_) -> void
    {
        close();

        options = options_;

        if(options.name.empty())
            options.name = "plot" + std::to_string(std::rand());

        if(options.legend.empty())
            options.legend = options.y;

        dataname = options.name + ".dat";
        plotname = options.name + ".plt";

        datafile.open(dataname);
        plotfile.open(plotname);

        // Output the header of each column
        datafile << std::left << std::setw(20) << options.x;
        for(auto y : options.y)
            datafile << std::left << std::setw(20) << y;
        datafile << std::endl;

        // Setup the Gnuplot script file
        plotfile << options.config << std::endl;

        // Define a list of titles for the curves
        plotfile << "titles = '";
        for(const auto& title : options.legend)
            plotfile << title << (&title != &options.legend.back() ? " " : "");
        plotfile << "'\n" << std::endl;

        // Define auxiliary variables for the plot
        auto imax = 1 + options.y.size();
        auto wait = 1.0/options.frequency;
        auto command = "gnuplot -persist " + plotname;

        // Write the lines for plotting the data
        plotfile << "pause " << wait << std::endl;
        plotfile << "plot for [i=2:" << imax << "] '" << dataname <<
            "' using 1:i with lines lt i-1 lw 2 title word(titles, i-1)\n" << std::endl;
        plotfile << "reread\n" << std::endl;

        // Open the Gnuplot plot
        pipe = popen(command.c_str(), "w");
    }

    auto close() -> void
    {
        if(pipe != nullptr)
            pclose(pipe);
        pipe = nullptr;
    }

    auto update(const ChemicalState& state, double t) -> void
    {
        quantity.update(state, t);
        datafile << std::left << std::setw(20) << quantity.value(options.x);
        for(auto y : options.y)
            datafile << std::left << std::setw(20) << quantity.value(y);
        datafile << std::endl;
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

auto ChemicalPlot::open(const ChemicalPlotOptions& options) -> void
{
    pimpl->open(options);
}

auto ChemicalPlot::update(const ChemicalState& state, double t) -> void
{
    pimpl->update(state, t);
}

} // namespace Reaktoro
