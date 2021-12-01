// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <memory>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>

namespace Reaktoro {

/// Use this class for storing the collection of states and their properties
class ChemicalField
{
public:
    using Iterator = std::vector<ChemicalState>::iterator;
    using ConstIterator = std::vector<ChemicalState>::const_iterator;

    /// Construct a ChemicalField instance when the number of the states(cells) and system is provided.
    ChemicalField(Index size, const ChemicalSystem& system);

    /// Construct a ChemicalField instance when the number of the states(cells) and an instance of one state.
    ChemicalField(Index size, const ChemicalState& state);

    /// Construct a copy of a ChemicalField instance.
    ChemicalField(const ChemicalField& other);

    /// Destroy this ChemicalField instance.
    virtual ~ChemicalField();

    /// Assign a copy of an ChemicalField instance.
    auto operator=(ChemicalField other) -> ChemicalField&;

    /// Return the number of states in chemical field.
    auto size() const -> Index;

    /// Return begin const iterator of this ChemicalField instance (for STL compatibility reasons).
    auto begin() const -> ConstIterator;

    /// Return begin iterator of this ChemicalField instance (for STL compatibility reasons).
    auto begin() -> Iterator;

    /// Return end const iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() const -> ConstIterator;

    /// Return end iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() -> Iterator;

    /// Return the ChemicalState object with given index.
    auto operator[](Index index) const -> const ChemicalState&;

    /// Return the ChemicalState object with given index.
    auto operator[](Index index) -> ChemicalState& ;

    /// Set all the states of the field by the one provided.
    auto set(const ChemicalState& state) -> void;

    /// Fetch the vector with temperatures in the states.
    auto temperature(ArrayXr values) -> void;

    /// Fetch the vector with pressures in the states.
    auto pressure(ArrayXr values) -> void;

    /// Fetch the vector with element amounts in the states.
    auto elementAmounts(ArrayXr values) -> void;

    // Output the quantities to the file
    auto output(const std::string& filename, const StringList& quantities) -> void;

    // Return the collection with properties
    auto properties() -> std::vector<ChemicalProps>&;

    // Return the collection with chemical state
    auto states() -> std::vector<ChemicalState>&;

private:

    struct Impl;

    std::unique_ptr<Impl> pimpl;

};

} // namespace Reaktoro
