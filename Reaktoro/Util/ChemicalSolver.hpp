// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.
//
//#pragma once
//
//// C++ includes
//#include <memory>
//#include <string>
//
//// Reaktoro includes
//#include <Reaktoro/Common/Index.hpp>
//#include <Reaktoro/Common/Matrix.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class ChemicalField;
//class ChemicalSystem;
//class ChemicalState;
//class Partition;
//class ReactionSystem;
//
///// A type that describes a solver for many chemical calculations.
//class ChemicalSolver
//{
//public:
//    // Forward declaration of the Array type.
//    template<typename T>
//    class Array;
//
//    // Forward declaration of the Grid type.
//    template<typename T>
//    class Grid;
//
//    /// Construct a default ChemicalSolver instance.
//    ChemicalSolver();
//
//    /// Construct a ChemicalSolver instance with given chemical system and field number of points.
//    ChemicalSolver(const ChemicalSystem& system, Index npoints);
//
//    /// Construct a ChemicalSolver instance with given reaction system and field number of points.
//    ChemicalSolver(const ReactionSystem& reactions, Index npoints);
//
//    /// Return the number of field points.
//    auto numPoints() const -> Index;
//
//    /// Return the number of equilibrium elements.
//    auto numEquilibriumElements() const -> Index;
//
//    /// Return the number of kinetic species.
//    auto numKineticSpecies() const -> Index;
//
//    /// Return the number of chemical components.
//    auto numComponents() const -> Index;
//
//    /// Set the partitioning of the chemical system.
//    auto setPartition(const Partition& partition) -> void;
//
//    /// Set the chemical state of all field points uniformly.
//    /// @param state The state of the chemical system.
//    auto setStates(const ChemicalState& state) -> void;
//
//    /// Set the chemical state of all field points.
//    /// @param states The array of states of the chemical system.
//    auto setStates(const Array<ChemicalState>& states) -> void;
//
//    /// Set the chemical state at a specified field point.
//    /// @param ipoint The index of the field point.
//    /// @param state The state of the chemical system.
//    auto setStateAt(Index ipoint, const ChemicalState& state) -> void;
//
//    /// Set the same chemical state at all specified field points.
//    /// @param ipoints The indices of the field points.
//    /// @param state The state of the chemical system.
//    auto setStateAt(const Array<Index>& ipoints, const ChemicalState& state) -> void;
//
//    /// Set the chemical state at all specified field points.
//    /// @param ipoints The indices of the field points.
//    /// @param states The states of the chemical system.
//    auto setStateAt(const Array<Index>& ipoints, const Array<ChemicalState>& states) -> void;
//
//    /// Equilibrate the chemical state at every field point.
//    auto equilibrate(Array<double> T, Array<double> P, Array<double> be) -> void;
//
//    /// Equilibrate the chemical state at every field point.
//    auto equilibrate(Array<double> T, Array<double> P, Grid<double> be) -> void;
//
//    /// React the chemical state at every field point.
//    auto react(double t, double dt) -> void;
//
//    /// Return the chemical state at given index.
//    auto state(Index i) const -> const ChemicalState&;
//
//    /// Return the chemical states at all field points.
//    auto states() const -> const std::vector<ChemicalState>&;
//
//    /// Return the molar amounts of the chemical components at every field point (in units of mol).
//    auto componentAmounts() -> const std::vector<VectorXr>&;
//
//    /// Return the molar amounts of each equilibrium species and their derivatives at every field point.
//    auto equilibriumSpeciesAmounts() -> const std::vector<ChemicalField>&;
//
//    /// Return the porosity at every field point.
//    auto porosity() -> const ChemicalField&;
//
//    /// Return the saturations of the fluid phases at every field point.
//    auto fluidSaturations() -> const std::vector<ChemicalField>&;
//
//    /// Return the densities of the fluid phases at every field point (in units of kg/m3).
//    auto fluidDensities() -> const std::vector<ChemicalField>&;
//
//    /// Return the volumes of the fluid phases at every field point (in units of m3).
//    auto fluidVolumes() -> const std::vector<ChemicalField>&;
//
//    /// Return the total volume of the fluid phases at every field point (in units of m3).
//    auto fluidTotalVolume() -> const ChemicalField&;
//
//    /// Return the total volume of the solid phases at every field point (in units of m3).
//    auto solidTotalVolume() -> const ChemicalField&;
//
//    /// Return the kinetic rates of the chemical components at every field point (in units of mol/s).
//    auto componentRates() -> const std::vector<ChemicalField>&;
//
//private:
//    struct Impl;
//
//    std::shared_ptr<Impl> pimpl;
//
//public:
//    /// A type that represents a one-dimensional array of values.
//    template<typename T>
//    class Array
//    {
//    public:
//        /// Construct a default Array instance.
//        Array()
//        {}
//
//        /// Construct a custom Array instance.
//        Array(const T* data, Index size)
//        : data(data), size(size) {}
//
//        /// Construct an Array instance from a vector-like instance.
//        /// The type of the vector must have public methods `data` and `size`.
//        template<typename VectorType>
//        Array(const VectorType& vec)
//        : data(vec.data()), size(vec.size()) {}
//
//        friend class ChemicalSolver;
//
//    private:
//        /// The pointer to a one-dimensional array of values.
//        const T* data = nullptr;
//
//        /// The size of the array.
//        Index size = 0;
//    };
//
//    /// A type that represents a two-dimensional array of values.
//    template<typename T>
//    class Grid
//    {
//    public:
//        /// Construct a default Grid instance.
//        Grid()
//        {}
//
//        /// Construct a custom Grid instance.
//        Grid(const T** data, Index rows, Index cols)
//        : data(data), rows(rows), cols(cols) {}
//
//        /// Construct an Grid instance from a vector-like instance.
//        /// The type of the vector must have public methods `data` and `size`.
//        template<typename VectorType>
//        Grid(const std::vector<VectorType>& vec)
//        {
//            pointers.reserve(vec.size());
//            for(const VectorType& v : vec)
//                pointers.push_back(v.data());
//            rows = vec.size();
//            cols = vec.front().size();
//        }
//
//        friend class ChemicalSolver;
//
//    private:
//        /// The pointer to a two-dimensional array of values.
//        const T** data = nullptr;
//
//        /// The pointer to a two-dimensional array of values.
//        std::vector<const T*> pointers;
//
//        /// The number of rows of the grid.
//        Index rows = 0;
//
//        /// The number of columns of the grid.
//        Index cols = 0;
//    };
//};
//
//} // namespace Reaktoro
//
