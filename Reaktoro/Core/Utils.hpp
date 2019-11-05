// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// Return the names of the entries in a container.
template<typename NamedValues>
auto names(const NamedValues& values) -> std::vector<std::string>;

/// Return the electrical charges of all species in a list of species
template<typename ChargedValues>
auto charges(const ChargedValues& values) -> Vector;

/// Return the molar masses of all species in a list of species (in units of kg/mol)
template<typename SpeciesValues>
auto molarMasses(const SpeciesValues& species) -> Vector;

/// Return the mole fractions of the species.
inline auto moleFractions(Composition n) -> ChemicalVector
{
    const auto nspecies = n.size();
    if(nspecies == 1)
        return ones(n);
    const ChemicalScalar nt = sum(n);
    if(nt != 0.0)
        return n / nt;
    else
        return zeros(n);
}

} // namespace Reaktoro

#include "Utils.hxx"
