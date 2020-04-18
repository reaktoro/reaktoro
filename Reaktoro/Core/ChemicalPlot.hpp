// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
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

// #pragma once

// // C++ includes
// #include <memory>
// #include <vector>
// #include <sstream>
// #include <string>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalState;
// class ChemicalSystem;
// class ReactionSystem;
// class StringList;

// /// A class used to create plots from sequence of chemical states.
// class ChemicalPlot
// {
// public:
//     /// Construct a default ChemicalPlot instance.
//     ChemicalPlot();

//     /// Construct a ChemicalPlot instance using a ChemicalSystem instance.
//     explicit ChemicalPlot(const ChemicalSystem& system);

//     /// Construct a ChemicalPlot instance using a ReactionSystem instance.
//     explicit ChemicalPlot(const ReactionSystem& reactions);

//     /// Destroy this ChemicalPlot instance.
//     virtual ~ChemicalPlot();

//     /// Set the name of the plot and the file names exported.
//     auto name(std::string name) -> void;

//     /// Set the quantity to be plotted along the x-axis.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.x("pH");
//     /// ~~~
//     /// @see ChemicalQuantity
//     auto x(std::string quantity) -> void;

//     /// Add a quantity to be plotted along the y-axis.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.y("elementAmount(Ca)");
//     /// plot.y("speciesMass(Calcite)");
//     /// ~~~
//     /// @note This method can be called multiple times.
//     /// @param quantity The quantity to be plotted.
//     /// @see ChemicalQuantity
//     auto y(std::string quantity) -> void;

//     /// Add a quantity to be plotted along the y-axis with a label.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.y("elementAmount(Ca)", "Ca [mol]");
//     /// plot.y("speciesMass(Calcite)", "Calcite [kg]");
//     /// ~~~
//     /// @note This method can be called multiple times.
//     /// @param quantity The quantity to be plotted.
//     /// @param label The label of the quantity displayed in the legend.
//     /// @see ChemicalQuantity
//     auto y(std::string quantity, std::string label) -> void;

//     /// Add discrete points in the plot.
//     /// **Usage Example**
//     /// @param label The label used in the legend to describe the points.
//     /// @param xpoints The x-coordinates of the points.
//     /// @param ypoints The y-coordinates of the points.
//     auto points(std::string label, std::vector<double> xpoints, std::vector<double> ypoints) -> void;

//     /// Add discrete points in the plot.
//     /// @param label The label used in the legend to describe the points.
//     /// @param xpoints The x-coordinates of the points separated by comma or space.
//     /// @param ypoints The y-coordinates of the points separated by comma or space.
//     auto points(std::string label, std::string xpoints, std::string ypoints) -> void;

//     /// Set the legend options.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.legend("left top");
//     /// plot.legend("right center");
//     /// ~~~
//     /// @see showlegend
//     auto legend(std::string) -> void;

//     /// Set `true` if legend should be displayed in the plot.
//     auto showlegend(bool active) -> void;

//     /// Return `true` if legend should be displayed in the plot.
//     auto showlegend() const -> bool;

//     /// Set the title of the plot.
//     auto title(std::string title) -> void;

//     /// Set the label of the x-axis.
//     /// @see ylabel
//     auto xlabel(std::string) -> void;

//     /// Set the label of the y-axis.
//     /// @see xlabel
//     auto ylabel(std::string) -> void;

//     /// Set the tics of the x-axis.
//     /// @see ytics
//     auto xtics(std::string) -> void;

//     /// Set the tics of the y-axis.
//     /// @see xtics
//     auto ytics(std::string) -> void;

//     /// Set the numeric display format of the x-axis.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.xformat("%f"); // sets floating point notation.
//     /// plot.xformat("%e"); // sets exponential notation using lower case `e`.
//     /// plot.xformat("%E"); // sets exponential notation using upper case `E`.
//     /// plot.xformat("%g"); // sets exponential notation like "%e", but shorter.
//     /// plot.xformat("%G"); // sets exponential notation like "%E", but shorter.
//     /// ~~~
//     /// @note The full list of supported format specifiers can be found in the [user's manual][gnuplot]
//     /// of Gnuplot v5.0, at page 123.
//     /// [gnuplot]: http://www.gnuplot.info/docs_5.0/gnuplot.pdf
//     /// @see yformat
//     auto xformat(std::string) -> void;

//     /// Set the numeric display format of the y-axis.
//     /// **Usage Example**
//     /// ~~~{.cpp}
//     /// plot.yformat("%f"); // sets floating point notation.
//     /// plot.yformat("%e"); // sets exponential notation using lower case `e`.
//     /// plot.yformat("%E"); // sets exponential notation using upper case `E`.
//     /// plot.yformat("%g"); // sets exponential notation like "%e", but shorter.
//     /// plot.yformat("%G"); // sets exponential notation like "%E", but shorter.
//     /// ~~~
//     /// @note The full list of supported format specifiers can be found in the [user's manual][gnuplot]
//     /// of Gnuplot v5.0, at page 123.
//     /// [gnuplot]: http://www.gnuplot.info/docs_5.0/gnuplot.pdf
//     /// @see xformat
//     auto yformat(std::string) -> void;

//     /// Set the x-axis to log-scale.
//     auto xlogscale(int base=10) -> void;

//     /// Set the y-axis to log-scale.
//     auto ylogscale(int base=10) -> void;

//     /// Set the refresh rate of the real-time plot.
//     auto frequency(unsigned frequency) -> void;

//     /// Inject a gnuplot command to the script file.
//     auto operator<<(std::string command) -> ChemicalPlot&;

//     /// Inject a gnuplot command to the script file.
//     auto operator<<(std::stringstream command) -> ChemicalPlot&;

//     /// Open the plot.
//     auto open() -> void;

//     /// Update the plot with a new chemical state and a tag.
//     auto update(const ChemicalState& state, double t) -> void;

//     /// Compare a ChemicalPlot instance for equality
//     auto operator==(const ChemicalPlot& other) -> bool;

// private:
//     struct Impl;

//     std::shared_ptr<Impl> pimpl;
// };

// } // namespace Reaktoro
