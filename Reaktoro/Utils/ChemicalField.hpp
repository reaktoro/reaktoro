// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Partition.hpp>

namespace Reaktoro {

// Forward declarations
class EquilibriumSensitivity;

/// A type that contains the values of a scalar field and its derivatives.
class ChemicalField
{
public:
    /// Construct a default ChemicalField instance.
    ChemicalField();

    /// Construct a ChemicalField instance with given chemical system partition.
    /// @param npoints The number of points in the field.
    /// @param partition The partition of the chemical system.
    ChemicalField(const Partition& partition, Index npoints);

    /// Set the field at the i-th point with a ChemicalScalar instance.
    /// @param i The index of the field point.
    /// @param scalar The chemical scalar to be set at the i-th point.
    /// @param sensitivity The equilibrium sensitivity at the i-th point.
    auto set(Index i, const ChemicalScalar& scalar, const EquilibriumSensitivity& sensitivity) -> void;

    /// Return size of the chemical field.
    auto size() const -> Index;

    /// The values of the scalar chemical field.
    Vector val;

    /// The derivatives of the scalar chemical field with respect to temperature.
    Vector ddT;

    /// The derivatives of the scalar chemical field with respect to pressure.
    Vector ddP;

    /// The derivatives of the scalar chemical field with respect to the amounts of each equilibrium element.
    std::vector<Vector> ddbe;

    /// The derivatives of the scalar chemical field with respect to the amounts of each kinetic species.
    std::vector<Vector> ddnk;

private:
    /// The partition of the chemical system.
    Partition partition;

    /// The number of points in the field.
    Index npoints = 0;

    /// Auxiliary vectors to avoid recurrent memory allocation.
    Vector scalar_ne, scalar_nk, scalar_be;
};

} // namespace Reaktoro

