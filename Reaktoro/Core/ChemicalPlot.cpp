// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright © 2014-2021 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "ChemicalPlot.hpp"

// // C++ includes
// #include <cstdio>
// #include <fstream>
// #include <iomanip>

// // Boost includes
// #include <boost/format.hpp>

// // Reaktoro includes
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Common/StringList.hpp>
// #include <Reaktoro/Common/StringUtils.hpp>
// #include <Reaktoro/Common/Units.hpp>
// #include <Reaktoro/Core/ChemicalQuantity.hpp>
// #include <Reaktoro/Core/ChemicalState.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>
// #include <Reaktoro/Core/ReactionSystem.hpp>

// // Ensure appropriate popen or pclose calls when compiling with MSVC
// #ifdef _MSC_VER
// #define popen _popen
// #define pclose _pclose
// #endif

// namespace Reaktoro {
// namespace {

// // Define the string that represents the preamble of the Gnuplot script
// const std::string gnuplot_preamble = R"(
// # Change the font
// set termoption enhanced
// set termoption font "Georgia,10"

// # Allow the use of macros
// set macros

// # Set a thick border
// set border lw 4

// # Set the ratio height/width to 0.618, the reciprocal of the golden ratio
// set size ratio 0.618

// # The line styles
// set style line 1 lt 2 lw 6 lc rgb '#0072bd' # blue
// set style line 2 lt 2 lw 6 lc rgb '#d95319' # orange
// set style line 3 lt 2 lw 6 lc rgb '#edb120' # yellow
// set style line 4 lt 2 lw 6 lc rgb '#7e2f8e' # purple
// set style line 5 lt 2 lw 6 lc rgb '#77ac30' # green
// set style line 6 lt 2 lw 6 lc rgb '#4dbeee' # light-blue
// set style line 7 lt 2 lw 6 lc rgb '#a2142f' # red

// # Allow a grid to be drawn
// set grid

// # Additional user options
// )";

// // Define the formatted string that represents the plot part of the Gnuplot script
// const std::string gnuplot_plot = R"xyz(
// # Defining the plot command
// COMMAND = "plot for [i=2:words(ytitles)+1] '%2%' using 1:i with lines ls i-1 title word(ytitles, i-1), \
//                 for [j=1:words(pfiles)] word(pfiles, j) using 1:2 with points ls j title word(ptitles, j)"

// # Check if the 'current' variable was defined as a gnuplot command-line parameter
// if(!exist('current')) @COMMAND; exit gnuplot

// # If 'current' is defined, then start the plotting loop
// previous = current
// current = system('%1% %2%')
// finished = system('%3%')
// pause %4%
// if(current ne previous && previous ne '') @COMMAND
// if(finished == 0) reread)xyz";

// } // namespace

// struct ChemicalPlot::Impl
// {
//     /// The chemical system instance
//     ChemicalSystem system;

//     /// The reaction system instance
//     ReactionSystem reactions;

//     /// The chemical quantity instance
//     ChemicalQuantity quantity;

//     /// The name of the plot and the name of the files output during the plot.
//     std::string name;

//     /// The quantity that spans the x-axis
//     std::string x = "t";

//     // The y data as pairs (legend, y-quantity)
//     std::vector<std::tuple<std::string, std::string>> y;

//     // The points data as triplets (legend, xpoints, ypoints).
//     std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> points;

//     /// The title of the plot.
//     std::string xlabel = "t";

//     /// The label of the y-axis.
//     std::string ylabel;

//     /// The boolean flag that indicates if the legend should be hidden
//     bool nolegend = false;

//     /// The Gnuplot commands used to configure the plot.
//     std::string config;

//     /// The frequency in which the plot is refreshed per second.
//     unsigned frequency = 30;

//     /// The name of the data file.
//     std::string dataname;

//     /// The name of the gnuplot script file.
//     std::string plotname;

//     /// The name of the file that is open to signal Gnuplot to start rereading.
//     std::string endname;

//     /// The output stream of the data file.
//     std::ofstream datafile;

//     /// The output stream of the gnuplot script file.
//     std::ofstream plotfile;

//     /// The file that is open to signal Gnuplot to start rereading.
//     std::ofstream endfile;

//     /// The pointer to the pipe connecting to Gnuplot
//     FILE* pipe = nullptr;

//     /// The iteration number for every update call
//     Index iteration = 0;

//     /// The counter of ChemicalPlot instances
//     static unsigned counter;

//     /// The ID of this ChemicalPlot instance (by order of creation)
//     unsigned id;

//     Impl()
//     : quantity(system)
//     {
//         id = counter++;
//     }

//     Impl(const ChemicalSystem& system)
//     : system(system), quantity(system)
//     {
//         id = counter++;
//     }

//     Impl(const ReactionSystem& reactions)
//     : system(reactions.system()), reactions(reactions), quantity(reactions)
//     {
//         id = counter++;
//     }

//     ~Impl()
//     {
//         close();
//     }

//     auto open() -> void
//     {
//         // Ensure the plot is closed
//         close();

//         // Make sure name is not empty
//         if(name.empty())
//             name = "plot" + std::to_string(id);

//         // Initialize the names of the data and gnuplot script files
//         dataname = name + ".dat";
//         plotname = name + ".plt";
//         endname  = name + ".end";

//         // Open the data and gnuplot script files
//         datafile.open(dataname);
//         plotfile.open(plotname);

//         // Set the output at higher precision to ensure smoother plots
//         datafile << std::setprecision(10);

//         // Output the name of each quantity in the data file
//         datafile << std::left << std::setw(20) << x;
//         for(auto item : y)
//             datafile << std::left << std::setw(20) << std::get<1>(item);
//         datafile << std::endl;

//         // For each discrete point data, output a file with given data
//         for(auto item : points)
//         {
//             // Auxiliary variables
//             const auto& legend  = std::get<0>(item);
//             const auto& xpoints = std::get<1>(item);
//             const auto& ypoints = std::get<2>(item);

//             // Ouput the data points to a file named `name-legend.dat`
//             std::ofstream file(name + "-" + legend + ".dat");

//             // Output the headings of the file
//             file << std::left << std::setw(20) << xlabel;
//             file << std::left << std::setw(20) << legend;
//             file << std::endl;

//             // Output each line of data (x, y)
//             for(Index i = 0; i < xpoints.size(); ++i)
//             {
//                 file << std::left << std::setw(20) << xpoints[i];
//                 file << std::left << std::setw(20) << ypoints[i];
//                 file << std::endl;
//             }
//         }

//         // Initialize the Gnuplot script file with the default preamble
//         plotfile << gnuplot_preamble;

//         // Apply the custom configuration
//         plotfile << config;

//         // Define Gnuplot a variable containing a list of titles for the lines plotted along the y-axis.
//         plotfile << "# The titles of each line plotted along the y-axis" << std::endl;
//         plotfile << "ytitles = \""; // e.g., ytitles = "legend1 legend2 legend3"
//         for(unsigned i = 0; i < y.size(); ++i)
//             plotfile << (i == 0 ? "" : " ") << "'" << std::get<0>(y[i]) << "'";
//         plotfile << "\"\n" << std::endl;

//         // Define Gnuplot a variable containing a list of titles for the discrete data points plotted along the y-axis.
//         plotfile << "# The titles of each discrete point data" << std::endl;
//         plotfile << "ptitles = \""; // e.g., ptitles = "legend1 legend2 legend3"
//         for(unsigned i = 0; i < points.size(); ++i)
//             plotfile << (i == 0 ? "" : " ") << "'" << std::get<0>(points[i]) << "'";
//         plotfile << "\"\n" << std::endl;

//         // Define Gnuplot a variable containing a list of file names containing the discrete data points.
//         plotfile << "# The names of the files containing discrete point data" << std::endl;
//         plotfile << "pfiles = \""; // e.g., pfiles = "plot-legend1.dat plot-legend2.dat plot-legend3.dat"
//         for(unsigned i = 0; i < points.size(); ++i)
//             plotfile << (i == 0 ? "" : " ") << name + "-" + std::get<0>(points[i]) + ".dat";
//         plotfile << "\"\n" << std::endl;

//         // On Windows, use the `dir` command on the data file to check its state.
//         // On any other OS, use the `ls -l` command instead.
// #if _WIN32
//         std::string file_status_cmd = "dir";
//         std::string file_exists_cmd = "if exist " + endname + " (echo 1) else (echo 0)";
// #else
//         std::string file_status_cmd = "ls -l";
//         std::string file_exists_cmd = "[ ! -e " + endname + " ]; echo $?";
// #endif
//         // Define auxiliary variables for the plot
//         auto wait = 1.0/frequency;

//         // Finalize the Gnuplot script
//         plotfile << boost::format(gnuplot_plot) % file_status_cmd % dataname % file_exists_cmd % wait;

//         // Flush the plot file to ensure its correct state before the plot starts
//         plotfile.flush();
//     }

//     auto close() -> void
//     {
//         if(pipe != nullptr)
//         {
//             // Create the file that signals Gnuplot to stop rereading the input script
//             endfile.open(endname);

//             // Close the pipe
//             pclose(pipe);

//             // Close the previously created file
//             endfile.close();

//             // Delete the end file
//             std::remove(endname.c_str());

//             // Set pipe to nullptr
//             pipe = nullptr;
//         }
//     }

//     auto update(const ChemicalState& state, double t) -> void
//     {
//         // Output the current chemical state to the data file.
//         quantity.update(state, t);
//         datafile << std::left << std::setw(20) << quantity.value(x);
//         for(auto item : y)
//         {
//             std::string qstr = std::get<1>(item);
//             auto val = (qstr == "i") ? iteration : quantity.value(qstr);
//             datafile << std::left << std::setw(20) << val;
//         }
//         datafile << std::endl;

//         // Open the Gnuplot plot after the first data has been output to the data file.
//         // This ensures that Gnuplot opens the plot without errors/warnings.
//         if(pipe == nullptr)
//         {
//             std::string command = ("gnuplot -persist -e \"current=''\" " + plotname + " >> gnuplot.log 2>&1");
//             pipe = popen(command.c_str(), "w");
//         }

//         // Update the iteration number
//         ++iteration;
//     }
// };

// // Initialize the counter of ChemicalPlot instances
// unsigned ChemicalPlot::Impl::counter = 0;

// ChemicalPlot::ChemicalPlot()
// : pimpl(new Impl())
// {}

// ChemicalPlot::ChemicalPlot(const ChemicalSystem& system)
// : pimpl(new Impl(system))
// {}

// ChemicalPlot::ChemicalPlot(const ReactionSystem& reactions)
// : pimpl(new Impl(reactions))
// {}

// ChemicalPlot::~ChemicalPlot()
// {}

// auto ChemicalPlot::name(std::string name) -> void
// {
//     pimpl->name = name;
// }

// auto ChemicalPlot::x(std::string quantity) -> void
// {
//     pimpl->x = quantity;
// }

// auto ChemicalPlot::y(std::string quantity) -> void
// {
//     pimpl->y.emplace_back(quantity, quantity);
// }

// auto ChemicalPlot::y(std::string quantity, std::string label) -> void
// {
//     pimpl->y.emplace_back(label, quantity);
// }

// auto ChemicalPlot::points(std::string label, std::vector<double> xpoints, std::vector<double> ypoints) -> void
// {
//     Assert(xpoints.size() == ypoints.size(), "Could not set the given points in the plot.",
//         "The number of x and y points do not match.");
//     pimpl->points.emplace_back(label, xpoints, ypoints);
// }

// auto ChemicalPlot::points(std::string label, std::string xpoints, std::string ypoints) -> void
// {
//     std::vector<double> xps, yps;
//     for(auto item : split(xpoints, ", "))
//         xps.push_back(tofloat(item));
//     for(auto item : split(ypoints, ", "))
//         yps.push_back(tofloat(item));
//     points(label, xps, yps);
// }

// auto ChemicalPlot::legend(std::string options) -> void
// {
//     *this << "set key " + options;
// }

// auto ChemicalPlot::showlegend(bool active) -> void
// {
//     pimpl->nolegend = !active;
//     *this << (active ? "set key on" : "set key off");
// }

// auto ChemicalPlot::showlegend() const -> bool
// {
//     return !pimpl->nolegend;
// }

// auto ChemicalPlot::title(std::string title) -> void
// {
//     *this << "set title '{/:Bold {" + title + "}}'";
// }

// auto ChemicalPlot::xlabel(std::string str) -> void
// {
//     pimpl->xlabel = str;
//     *this << "set xlabel '" + str + "'";
// }

// auto ChemicalPlot::ylabel(std::string str) -> void
// {
//     pimpl->ylabel = str;
//     *this << "set ylabel '" + str + "'";
// }

// auto ChemicalPlot::xtics(std::string str) -> void
// {
//     *this << "set xtics " + str;
// }

// auto ChemicalPlot::ytics(std::string str) -> void
// {
//     *this << "set ytics " + str;
// }

// auto ChemicalPlot::xformat(std::string str) -> void
// {
//     *this << "set format x '" + str + "'";
// }

// auto ChemicalPlot::yformat(std::string str) -> void
// {
//     *this << "set format y '" + str + "'";
// }

// auto ChemicalPlot::xlogscale(int base) -> void
// {
//     *this << "set logscale x " + std::to_string(base);
// }

// auto ChemicalPlot::ylogscale(int base) -> void
// {
//     *this << "set logscale y " + std::to_string(base);
// }

// auto ChemicalPlot::frequency(unsigned frequency) -> void
// {
//     pimpl->frequency = frequency;
// }

// auto ChemicalPlot::operator<<(std::string command) -> ChemicalPlot&
// {
//     pimpl->config.append(command + "\n");
//     return *this;
// }

// auto ChemicalPlot::operator<<(std::stringstream command) -> ChemicalPlot&
// {
//     pimpl->config.append(command.str());
//     return *this;
// }

// auto ChemicalPlot::open() -> void
// {
//     pimpl->open();
// }

// auto ChemicalPlot::update(const ChemicalState& state, double t) -> void
// {
//     pimpl->update(state, t);
// }

// auto ChemicalPlot::operator==(const ChemicalPlot& other) -> bool
// {
//     return pimpl == other.pimpl;
// }

// } // namespace Reaktoro
