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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>

namespace Reaktoro {

class ChemicalField
{
public:
    using Iterator = std::vector<ChemicalState>::iterator;

    using ConstIterator = std::vector<ChemicalState>::const_iterator;

    ChemicalField(Index size, const ChemicalSystem& system);

    ChemicalField(Index size, const ChemicalState& state);

    auto size() const -> Index { return m_size; }

    auto begin() const -> ConstIterator { return m_states.cbegin(); }

    auto begin() -> Iterator { return m_states.begin(); }

    auto end() const -> ConstIterator { return m_states.cend(); }

    auto end() -> Iterator { return m_states.end(); }

    auto operator[](Index index) const -> const ChemicalState& { return m_states[index]; }

    auto operator[](Index index) -> ChemicalState& { return m_states[index]; }

    auto set(const ChemicalState& state) -> void;

    auto temperature(VectorRef values) -> void;

    auto pressure(VectorRef values) -> void;

    auto elementAmounts(VectorRef values) -> void;

    auto output(std::string filename, StringList quantities) -> void;

private:
    /// The number of degrees of freedom in the chemical field.
    Index m_size;

//    Vector temperatures;
//
//    Vector pressures;
//
//    /// The matrix of amounts for every element (
//    Matrix element_amounts;

    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem m_system;

    /// The chemical states in the chemical field
    std::vector<ChemicalState> m_states;

    /// The chemical states in the chemical field
    std::vector<ChemicalProperties> m_properties;
};

} // namespace Reaktoro
