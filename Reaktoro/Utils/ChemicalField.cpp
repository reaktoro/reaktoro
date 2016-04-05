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

#include "ChemicalField.hpp"

// Reaktoro includes
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>

namespace Reaktoro {

ChemicalField::ChemicalField()
{}

ChemicalField::ChemicalField(const Partition& partition, Index npoints)
: val(npoints), ddT(npoints), ddP(npoints),
  ddbe(partition.numEquilibriumElements(), Vector(npoints)),
  ddnk(partition.numKineticSpecies(), Vector(npoints)),
  partition(partition), npoints(npoints)
{}

auto ChemicalField::set(Index i, const ChemicalScalar& scalar, const EquilibriumSensitivity& sensitivity) -> void
{
    // The indices of the equilibrium and kinetic species
    const Indices& ispecies_e = partition.indicesEquilibriumSpecies();
    const Indices& ispecies_k = partition.indicesKineticSpecies();

    // Auxiliary references to sensitivity values
    const Vector& ne_T  = sensitivity.dnedT;
    const Vector& ne_P  = sensitivity.dnedP;
    const Matrix& ne_be = sensitivity.dnedbe;

    // Extract the derivatives of scalar w.r.t. amounts of equilibrium species
    scalar_ne = rows(scalar.ddn, ispecies_e);

    // Extract the derivatives of scalar w.r.t. amounts of kinetic species
    scalar_nk = rows(scalar.ddn, ispecies_k);

    // Calculte the derivatives of scalar w.r.t. amounts of equilibrium elements
    scalar_be = tr(scalar_ne) * ne_be;

    // Set the i-th position of the scalar field with given scalar value
    val[i] = scalar.val;

    // Set derivative w.r.t. temperature at the i-th position
    ddT[i] = scalar.ddt + dot(scalar_ne, ne_T);

    // Set derivative w.r.t. pressure at the i-th position
    ddP[i] = scalar.ddp + dot(scalar_ne, ne_P);

    // Set derivative w.r.t. amounts of equilibrium elements at the i-th position
    for(Index j = 0; j < ddbe.size(); ++j)
        ddbe[j][i] = scalar_be[j];

    // Set derivative w.r.t. amounts of kinetic species at the i-th position
    for(Index j = 0; j < ddnk.size(); ++j)
        ddnk[j][i] = scalar_nk[j];
}

auto ChemicalField::size() const -> Index
{
    return npoints;
}

}  // namespace Reaktoro
