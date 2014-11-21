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
//#include "PartitionUtils.hpp"
//
//// Reaktor includes
//#include <Reaktor/Core/Multiphase.hpp>
//#include <Reaktor/Core/MultiphaseUtils.hpp>
//#include <Reaktor/Core/Partition.hpp>
//#include <Reaktor/Core/PartitionUtils.hpp>
//
//namespace Reaktor {
//
//auto numSpecies(const Partition& partition) -> unsigned
//{
//    return numEquilibriumSpecies(partition) +
//           numKineticSpecies(partition) +
//           numInertSpecies(partition);
//}
//
//auto numEquilibriumSpecies(const Partition& partition) -> unsigned
//{
//    return partition.equilibriumSpeciesIndices().size();
//}
//
//auto numKineticSpecies(const Partition& partition) -> unsigned
//{
//    return partition.kineticSpeciesIndices().size();
//}
//
//auto numInertSpecies(const Partition& partition) -> unsigned
//{
//    return partition.inertSpeciesIndices().size();
//}
//
//auto elementIndicesInEquilibriumSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return elementIndicesInSpecies(multiphase, partition.equilibriumSpeciesIndices());
//}
//
//auto elementIndicesInKineticSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return elementIndicesInSpecies(multiphase, partition.kineticSpeciesIndices());
//}
//
//auto elementIndicesInInertSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return elementIndicesInSpecies(multiphase, partition.inertSpeciesIndices());
//}
//
//auto phaseIndicesWithEquilibriumSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return phaseIndicesWithSpecies(multiphase, partition.equilibriumSpeciesIndices());
//}
//
//auto phaseIndicesWithKineticSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return phaseIndicesWithSpecies(multiphase, partition.kineticSpeciesIndices());
//}
//
//auto phaseIndicesWithInertSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
//{
//    return phaseIndicesWithSpecies(multiphase, partition.inertSpeciesIndices());
//}
//
//auto equilibriumRows(const Partition& partition, const Vector& vec) -> Vector
//{
//    return vec.elem(arma::uvec(partition.equilibriumSpeciesIndices()));
//}
//
//auto kineticRows(const Partition& partition, const Vector& vec) -> Vector
//{
//    return vec.elem(arma::uvec(partition.kineticSpeciesIndices()));
//}
//
//auto inertRows(const Partition& partition, const Vector& vec) -> Vector
//{
//    return vec.elem(arma::uvec(partition.inertSpeciesIndices()));
//}
//
//auto equilibriumCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//    return mat.cols(arma::uvec(partition.equilibriumSpeciesIndices()));
//}
//
//auto kineticCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//    return mat.cols(arma::uvec(partition.kineticSpeciesIndices()));
//}
//
//auto inertCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//    return mat.cols(arma::uvec(partition.inertSpeciesIndices()));
//}
//
//auto equilibriumRowsCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//	const arma::uvec indices = partition.equilibriumSpeciesIndices();
//    return mat.submat(indices, indices);
//}
//
//auto kineticRowsCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//	const arma::uvec indices = partition.kineticSpeciesIndices();
//    return mat.submat(indices, indices);
//}
//
//auto inertRowsCols(const Partition& partition, const Matrix& mat) -> Matrix
//{
//	const arma::uvec indices = partition.inertSpeciesIndices();
//    return mat.submat(indices, indices);
//}
//
//auto equilibriumFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix
//{
//    const arma::uvec& ispecies = partition.equilibriumSpeciesIndices();
//    const arma::uvec& ielements = elementIndicesInEquilibriumSpecies(multiphase, partition);
//    return mat.submat(ielements, ispecies);
//}
//
//auto kineticFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix
//{
//    const arma::uvec& ispecies = partition.kineticSpeciesIndices();
//    const arma::uvec& ielements = elementIndicesInKineticSpecies(multiphase, partition);
//    return mat.submat(ielements, ispecies);
//}
//
//auto inertFormulaMatrix(const Multiphase& multiphase, const Partition& partition, const Matrix& mat) -> Matrix
//{
//    const arma::uvec& ispecies = partition.inertSpeciesIndices();
//    const arma::uvec& ielements = elementIndicesInInertSpecies(multiphase, partition);
//    return mat.submat(ielements, ispecies);
//}
//
//} // namespace Reaktor
