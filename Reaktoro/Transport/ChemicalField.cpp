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

#include "ChemicalField.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

ChemicalField::ChemicalField(Index size, const ChemicalSystem& system)
: m_size(size),
  m_system(system),
  m_states(size, ChemicalState(system)),
  m_properties(size, ChemicalProperties(system))
{}

ChemicalField::ChemicalField(Index size, const ChemicalState& state)
: m_size(size),
  m_system(state.system()),
  m_states(size, state),
  m_properties(size, state.properties())
{}

auto ChemicalField::set(const ChemicalState& state) -> void
{
    for(auto& item : m_states)
        item = state;
}

auto ChemicalField::temperature(VectorRef values) -> void
{
    const Index len = size();
    for(Index i = 0; i < len; ++i)
        values[i] = m_states[i].temperature();
}

auto ChemicalField::pressure(VectorRef values) -> void
{
    const Index len = size();
    for(Index i = 0; i < len; ++i)
        values[i] = m_states[i].pressure();
}

auto ChemicalField::elementAmounts(VectorRef values) -> void
{
    const Index len = size();
    const Index num_elements = m_system.numElements();
    Index offset = 0;
    for(Index i = 0; i < len; ++i, offset += num_elements)
        values.segment(offset, num_elements) = m_states[i].elementAmounts();
}

auto ChemicalField::output(std::string filename, StringList quantities) -> void
{
    ChemicalOutput out(m_system);
    out.filename(filename);
    for(auto quantity : quantities)
        out.add(quantity);
}

} // namespace Reaktoro
