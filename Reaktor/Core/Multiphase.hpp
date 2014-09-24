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

#pragma once

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <memory>

namespace Reaktor {

// Forward declarations
class Phase;
class Species;

/// Provide a computational representation of a multiphase system
/// @see Species, Phase
/// @ingroup Core
class Multiphase
{
public:
    /// Construct a default ChemicalSystem instance
    Multiphase();

    /// Construct a ChemicalSystem instance with a given list of phases
    /// @param phases The list of Phase instances
    /// @see Phase
    explicit Multiphase(const std::vector<Phase>& phases);

    /// Get the names of the elements in the multiphase system
    auto elements() const -> const std::vector<std::string>&;

    /// Get the names of the species in the multiphase system
    auto species() const -> const std::vector<Species>&;

    /// Get the phases in the multiphase system
    auto phases() const -> const std::vector<Phase>&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a Phases instance
auto operator<<(std::ostream& out, const Multiphase& multiphase) -> std::ostream&;

} // namespace Reaktor
