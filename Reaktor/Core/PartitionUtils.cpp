// Reaktor is a C++ library for computational reaction modelling.
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

#include "PartitionUtils.hpp"

// Reaktor includes
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/MultiphaseUtils.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Core/PartitionUtils.hpp>

namespace Reaktor {

auto numSpecies(const Partition& partition) -> unsigned
{
    return numEquilibriumSpecies(partition) +
        numKineticSpecies(partition) +
        numInertSpecies(partition);
}

auto numEquilibriumSpecies(const Partition& partition) -> unsigned
{
    return partition.equilibriumSpeciesIndices().size();
}

auto numKineticSpecies(const Partition& partition) -> unsigned
{
    return partition.kineticSpeciesIndices().size();
}

auto numInertSpecies(const Partition& partition) -> unsigned
{
    return partition.inertSpeciesIndices().size();
}

auto indicesPhasesWithEquilibriumSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
{
    return phaseIndicesWithSpecies(multiphase, partition.equilibriumSpeciesIndices());
}

auto indicesPhasesWithKineticSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
{
    return phaseIndicesWithSpecies(multiphase, partition.kineticSpeciesIndices());
}

auto indicesPhasesWithInertSpecies(const Multiphase& multiphase, const Partition& partition) -> Indices
{
    return phaseIndicesWithSpecies(multiphase, partition.inertSpeciesIndices());
}

auto equilibriumRows(const Partition& partition, const Vector& vec) -> SubVector
{
    return vec.elem(arma::uvec(partition.equilibriumSpeciesIndices()));
}

auto kineticRows(const Partition& partition, const Vector& vec) -> SubVector
{
    return vec.elem(arma::uvec(partition.kineticSpeciesIndices()));
}

auto inertRows(const Partition& partition, const Vector& vec) -> SubVector
{
    return vec.elem(arma::uvec(partition.inertSpeciesIndices()));
}

auto equilibriumCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    return mat.cols(arma::uvec(partition.equilibriumSpeciesIndices()));
}

auto kineticCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    return mat.cols(arma::uvec(partition.kineticSpeciesIndices()));
}

auto inertCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    return mat.cols(arma::uvec(partition.inertSpeciesIndices()));
}

auto equilibriumRowsCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
	const arma::uvec indices = partition.equilibriumSpeciesIndices();
    return mat.submat(indices, indices);
}

auto kineticRowsCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
	const arma::uvec indices = partition.kineticSpeciesIndices();
    return mat.submat(indices, indices);
}

auto inertRowsCols(const Partition& partition, const Matrix& mat) -> SubMatrix
{
	const arma::uvec indices = partition.inertSpeciesIndices();
    return mat.submat(indices, indices);
}

auto equilibriumFormulaMatrix(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    const arma::uvec& ispecies = partition.equilibriumSpeciesIndices();
    const arma::uvec& ielements = partition.equilibriumElementIndices();
    return mat.submat(ielements, ispecies);
}

auto kineticFormulaMatrix(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    const arma::uvec& ispecies = partition.kineticSpeciesIndices();
    const arma::uvec& ielements = partition.kineticElementIndices();
    return mat.submat(ielements, ispecies);
}

auto inertFormulaMatrix(const Partition& partition, const Matrix& mat) -> SubMatrix
{
    const arma::uvec& ispecies = partition.inertSpeciesIndices();
    const arma::uvec& ielements = partition.inertElementIndices();
    return mat.submat(ielements, ispecies);
}

} // namespace Reaktor
