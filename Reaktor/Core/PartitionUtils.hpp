//// Reaktor is a C++ library for computational reaction modelling.
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
//// Reaktor includes
//#include <Reaktor/Common/Index.hpp>
//#include <Reaktor/Common/Vector.hpp>
//#include <Reaktor/Common/Matrix.hpp>
//
//namespace Reaktor {
//
//// Forward declarations
//class Multiphase;
//class Partition;
//
///// Get the number of species in a partition of a multiphase system
///// @param partition The partition of the multiphase system
//auto numSpecies(const Partition& partition) -> unsigned;
//
///// Get the number of equilibrium species in a partition of a multiphase system
///// @param partition The partition of the multiphase system
//auto numEquilibriumSpecies(const Partition& partition) -> unsigned;
//
///// Get the number of kinetic species in a partition of a multiphase system
///// @param partition The partition of the multiphase system
//auto numKineticSpecies(const Partition& partition) -> unsigned;
//
///// Get the number of inert species in a partition of a multiphase system
///// @param partition The partition of the multiphase system
//auto numInertSpecies(const Partition& partition) -> unsigned;
//
///// Get the indices of the elements in a multiphase system that compose the equilibrium species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto elementIndicesInEquilibriumSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the indices of the elements in a multiphase system that compose the kinetic species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto elementIndicesInKineticSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the indices of the elements in a multiphase system that compose the inert species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto elementIndicesInInertSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the indices of the phases in a multiphase system that contains equilibrium species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto phaseIndicesWithEquilibriumSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the indices of the phases in a multiphase system that contains kinetic species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto phaseIndicesWithKineticSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the indices of the phases in a multiphase system that contains inert species
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
//auto phaseIndicesWithInertSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices;
//
///// Get the rows of a vector that correspond to the equilibrium species
///// @param partition The partition of the multiphase system
///// @param vec The vector with values for all species
//auto equilibriumRows(const Partition& partition, const Vector& vec) -> Vector;
//
///// Get the rows of a vector that correspond to the kinetic species
///// @param partition The partition of the multiphase system
///// @param vec The vector with values for all species
//auto kineticRows(const Partition& partition, const Vector& vec) -> Vector;
//
///// Get the rows of a vector that correspond to the inert species
///// @param partition The partition of the multiphase system
///// @param vec The vector with values for all species
//auto inertRows(const Partition& partition, const Vector& vec) -> Vector;
//
///// Get the columns of a matrix that correspond to the equilibrium species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto equilibriumCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the columns of a matrix that correspond to the kinetic species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto kineticCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the columns of a matrix that correspond to the inert species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto inertCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of a matrix that correspond to the equilibrium species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto equilibriumRowsCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of a matrix that correspond to the kinetic species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto kineticRowsCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of a matrix that correspond to the inert species
///// @param partition The partition of the multiphase system
///// @param mat The matrix with values for all species
//auto inertRowsCols(const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of the formula matrix that corresponds to the elements and species in the equilibrium partition
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
///// @param mat The formula matrix of the multiphase system
//auto equilibriumFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of the formula matrix that corresponds to the elements and species in the kinetic partition
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
///// @param mat The formula matrix of the multiphase system
//auto kineticFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix;
//
///// Get the rows and columns of the formula matrix that corresponds to the elements and species in the inert partition
///// @param multiphase The multiphase system
///// @param partition The partition of the multiphase system
///// @param mat The formula matrix of the multiphase system
//auto inertFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix;
//
//} // namespace Reaktor
