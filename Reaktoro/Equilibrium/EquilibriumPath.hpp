// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2021 Allan Leal
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
// #include <string>
// #include <vector>

// // Reaktoro includes
// #include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
// #include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
// #include <Reaktoro/Math/ODE.hpp>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalOutput;
// class ChemicalPlot;
// class ChemicalSystem;
// class ChemicalState;
// class Partition;

// /// A struct that describes the options from an equilibrium path calculation.
// struct EquilibriumPathOptions
// {
//     /// The options for the chemical equilibrium calculations.
//     EquilibriumOptions equilibrium;

//     /// The options for the ODE solver
//     ODEOptions ode;

//     /// The maximum step length during the equilibrium path calculation.
//     double maxstep = 0.1;
// };

// /// A struct that describes the result of an equilibrium path calculation.
// struct EquilibriumPathResult
// {
//     /// The accumulated result of the equilibrium calculations.
//     EquilibriumResult equilibrium;
// };

// /// A class that describes a path of equilibrium states.
// class EquilibriumPath
// {
// public:
//     /// Construct an EquilibriumPath instance
//     explicit EquilibriumPath(const ChemicalSystem& system);

//     /// Construct a copy of an EquilibriumPath instance
//     EquilibriumPath(const EquilibriumPath& other);

//     /// Destroy an EquilibriumPath instance
//     virtual ~EquilibriumPath();

//     /// Assign an EquilibriumPath instance to this instance
//     auto operator=(EquilibriumPath other) -> EquilibriumPath&;

//     /// Set the options for the equilibrium path calculation and visualization
//     auto setOptions(const EquilibriumPathOptions& options) -> void;

//     /// Set the partition of the chemical system
//     auto setPartition(const Partition& partition) -> void;

//     /// Solve the path of equilibrium states between two chemical states
//     auto solve(const ChemicalState& state_i, const ChemicalState& state_f) -> EquilibriumPathResult;

//     /// Return a ChemicalPlot instance.
//     /// The returned ChemicalOutput instance must be properly configured
//     /// before the method EquilibriumPath::solve is called.
//     /// Changes in this ChemicalOutput instance are observed by the
//     /// EquilibriumPath object.
//     auto output() -> ChemicalOutput;

//     /// Return a ChemicalPlot instance.
//     /// The returned ChemicalPlot instance must be properly configured
//     /// before the method EquilibriumPath::solve is called.
//     /// Changes in this ChemicalPlot instance are observed by the
//     /// EquilibriumPath object.
//     auto plot() -> ChemicalPlot;

//     /// Return a collection of ChemicalPlot instances.
//     /// The returned ChemicalPlot instances must be properly configured
//     /// before the method EquilibriumPath::solve is called.
//     /// Changes in theses ChemicalPlot instances are observed by the
//     /// EquilibriumPath object.
//     auto plots(unsigned num) -> std::vector<ChemicalPlot>;

//     /// Return the chemical system in the equilibrium path definition.
//     auto system() const -> const ChemicalSystem&;

//     /// Return the partition of the chemical system in the equilibrium path definition.
//     auto partition() const -> const Partition&;

// private:
//     struct Impl;

//     std::unique_ptr<Impl> pimpl;
// };

// } // namespace Reaktoro
