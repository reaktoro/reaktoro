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

#include "ChemicalSystem.hpp"

// Reaktor includes
#include <Reaktor/Math/MathUtils.hpp>

namespace Reaktor {

struct ChemicalSystem::Impl
{
    /// The data used to construct the chemical system
    ChemicalSystemData data;

    /// The list of species in the chemical system
    SpeciesList species;

    /// The list of elements in the chemical system
    ElementList elements;

    Impl()
    {}

    Impl(const ChemicalSystemData& data)
    : data(data),
      species(collectSpecies(data.phases)),
      elements(collectElements(species))
    {}
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const ChemicalSystemData& data)
: pimpl(new Impl(data))
{}

auto ChemicalSystem::elements() const -> const ElementList&
{
    return pimpl->elements;
}

auto ChemicalSystem::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto ChemicalSystem::phases() const -> const PhaseList&
{
    return pimpl->data.phases;
}

auto ChemicalSystem::models() const -> const ChemicalSystemModels&
{
    return pimpl->data.models;
}

auto formulaMatrix(const ChemicalSystem& system) -> Matrix
{
    return formulaMatrix(system.species(), system.elements());
}

auto balanceMatrix(const ChemicalSystem& system) -> Matrix
{
    Matrix balance_matrix = formulaMatrix(system);
    Vector charges = collectCharges(system.species());
    balance_matrix = arma::join_vert(balance_matrix, charges);
    Indices components = linearlyIndependentRows(balance_matrix, balance_matrix);
    return balance_matrix;
}

} // namespace Reaktor
