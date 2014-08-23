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
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Core/Functions.hpp>

namespace Reaktor {

// Forward declarations
class Species;

/// Provide a computational representation of a phase
///
/// The Phase class can be seen as a collection of Species instances with
/// functionalities to compute their standard chemical potentials and activities.
///
/// Instances of class Phase are the building blocks of a Phases instance,
/// which computationaly represents multiphase system.
///
/// A Phase instance can be produced from an AqueousPhase, GaseousPhase, or
/// MineralPhase instance. These classes should be used together with the Database
/// class to define a phase -- specifying their species, their standard chemical
/// potential models, and their activity models.
///
/// @see Species, Phases
/// @ingroup Core
class Phase
{
public:
    /// Construct a default Phase instance
    Phase();

    /// Construct a copy of a Phase instance
    Phase(const Phase& other);

    /// Destroy this Phase instance
    virtual ~Phase();

    /// Assign a Phase instance to this instance
    auto operator=(Phase other) -> Phase&;

    /// Set the name of the phase
    /// @param name The name of the phase
    /// @return A reference to this Phase instance
    auto setName(const std::string& name) -> Phase&;

    /// Set the chemical species that compose the phase
    /// @param species The species that compose the phase
    /// @return A reference to this Phase instance
    /// @see Species
    auto setSpecies(const std::vector<Species>& species) -> Phase&;

    /// Set the concentration function of the phase
    /// @param concentration The concentration function of the phase
    /// @return A reference to this Phase instance
    /// @see Concentration
    auto setConcentration(const Concentration& concentration) -> Phase&;

    /// Get the name of the phase
    auto name() const -> const std::string&;

    /// Get the chemical species of the phase
    auto species() const -> const std::vector<Species>&;

    /// Get the concentration function of the phase
    auto concentration() const -> const Concentration&;

    /// Checks if this Phase instance is equal to another
    auto operator==(const Phase& phase) const -> bool;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output a Phase instance
auto operator<<(std::ostream& out, const Phase& phase) -> std::ostream&;

} /* namespace Reaktor */
