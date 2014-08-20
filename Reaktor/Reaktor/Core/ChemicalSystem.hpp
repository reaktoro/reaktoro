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

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class Phase;
class Species;

/// Provide a computational representation of a multiphase chemical system
///
/// The ChemicalSystem class is used to computationaly model a multiphase chemical system.
/// It is instantiated from a list of Phase instances defining the phases of the system.
/// Once it is built, it cannot be modified (i.e., it is an immutable object).
///
/// @see Phase, Species
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a ChemicalSystem instances
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with a given list of phases
    /// @param phases The list of phases that define the chemical system
    /// @see Phase
    explicit ChemicalSystem(const std::vector<Phase>& phases);

    /// Get the chemical species in the system
    auto species() const -> const std::vector<Species>&;

    /// Get the chemical species with given index
    /// @param i The index of the species
    auto species(const Index& i) const -> const Species&;

    /// Get the phases in the system
    auto phases() const -> const std::vector<Phase>&;

    /// Get the phase with given index
    /// @param i The index of the phase
    auto phase(const Index& i) const -> const Phase&;

    /// Get the chemical elements in the system
    auto elements() const -> const std::vector<std::string>&;

    /// Get the name of a chemical element in the system
    /// @param i The index of the chemical element
    auto element(const Index& i) const -> const std::string&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a ChemicalSystem instance
auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&;

} /* namespace Reaktor */
