//// Reaktoro is a C++ library for computational reaction modelling.
////
//// Copyright (C) 2014 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#pragma once
//
//// C++ includes
//#include <memory>
//#include <vector>
//#include <sstream>
//#include <string>
//
//namespace Reaktoro {
//
//// Forward declarations
//class ChemicalState;
//class ChemicalSystem;
//class ReactionSystem;
//
//class ChemicalPlot
//{
//public:
//    ChemicalPlot();
//
//    explicit ChemicalPlot(const ChemicalSystem& system);
//
//    explicit ChemicalPlot(const ReactionSystem& reactions);
//
//    virtual ~ChemicalPlot();
//
//    auto operator<<(std::string str) -> ChemicalPlot&;
//
//    auto operator<<(const std::stringstream& ss) -> ChemicalPlot&;
//
//    auto update(double t, const ChemicalState& state) -> void;
//
//    /// The strings representing the what to be plotted.
//    std::vector<std::string> what;
//
//    /// The strings representing the names of each curve.
//    std::vector<std::string> legend;
//
//    /// The integer that control how often the plot is refreshed.
//    unsigned every = 10;
//
//    /// The name of the file where the data is output
//    std::string filename;
//
//private:
//    struct Impl;
//
//    std::shared_ptr<Impl> pimpl;
//};
//
//} // namespace Reaktoro
