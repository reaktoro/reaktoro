//// Reaktoro is a unified framework for modeling chemically reactive systems.
////
//// Copyright (C) 2014-2015 Allan Leal
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
//// Reaktoro includes
//#include <Reaktoro/Common/ChemicalScalar.hpp>
//#include <Reaktoro/Math/Matrix.hpp>
//#include <Reaktoro/Core/Partition.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//struct EquilibriumSensitivity;
//
///// A type that contains the values of a scalar field and its derivatives.
//class ChemicalField
//{
//public:
//    /// Construct a default ChemicalField instance.
//    ChemicalField();
//
//    /// Construct a ChemicalField instance with given chemical system partition.
//    /// @param npoints The number of points in the field.
//    /// @param partition The partition of the chemical system.
//    ChemicalField(const Partition& partition, Index npoints);
//
//    /// Construct a copy of a ChemicalField instance.
//    ChemicalField(const ChemicalField& other);
//
//    /// Destroy this instance.
//    virtual ~ChemicalField();
//
//    /// Construct a copy of a ChemicalField instance.
//    auto operator=(ChemicalField other) -> ChemicalField&;
//
//    /// Set the field at the i-th point with a ChemicalScalar instance.
//    /// @param i The index of the field point.
//    /// @param scalar The chemical scalar to be set at the i-th point.
//    /// @param sensitivity The equilibrium sensitivity at the i-th point.
//    auto set(Index i, const ChemicalScalar& scalar, const EquilibriumSensitivity& sensitivity) -> void;
//
//    /// Return the partition of the chemical system.
//    auto partition() const -> const Partition&;
//
//    /// Return the size of the chemical field.
//    auto size() const -> Index;
//
//    /// Return a reference to the values of the chemical field.
//    auto val() -> VectorRef;
//
//    /// Return a const reference to the values of the chemical field.
//    auto val() const -> VectorConstRef;
//
//    /// Return a reference to the derivatives w.r.t. temperature of the chemical field.
//    auto ddT() -> VectorRef;
//
//    /// Return a const-reference to the derivatives w.r.t. temperature of the chemical field.
//    auto ddT() const -> VectorConstRef;
//
//    /// Return a reference to the derivatives w.r.t. pressure of the chemical field.
//    auto ddP() -> VectorRef;
//
//    /// Return a const-reference to the derivatives w.r.t. pressure of the chemical field.
//    auto ddP() const -> VectorConstRef;
//
//    /// Return a reference to the derivatives w.r.t. molar amounts of equilibrium elements of the chemical field.
//    auto ddbe() -> std::vector<Vector>&;
//
//    /// Return a const-reference to the derivatives w.r.t. molar amounts of equilibrium elements of the chemical field.
//    auto ddbe() const -> const std::vector<Vector>&;
//
//    /// Return a reference to the derivatives w.r.t. molar amounts of kinetic species of the chemical field.
//    auto ddnk() -> std::vector<Vector>&;
//
//    /// Return a const-reference to the derivatives w.r.t. molar amounts of kinetic species of the chemical field.
//    auto ddnk() const -> const std::vector<Vector>&;
//
//private:
//    struct Impl;
//
//    std::unique_ptr<Impl> pimpl;
//};
//
///// Output a ChemicalField instance.
//auto operator<<(std::ostream& out, const ChemicalField& f) -> std::ostream&;
//
//} // namespace Reaktoro
//
