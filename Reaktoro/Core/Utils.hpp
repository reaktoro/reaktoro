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
template<typename Derived>
auto molarFractions(const Eigen::MatrixBase<Derived>& n) -> ChemicalVector;

} // namespace Reaktoro

#include "Utils.hxx"
