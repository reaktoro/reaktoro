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

struct ChemicalField::Impl
{
    /// The number of degrees of freedom in the chemical field
    Index m_size;

    /// The chemical system common to all degrees of freedom in the chemical field
    ChemicalSystem m_system;

    /// The chemical states in the chemical field
    std::vector<ChemicalState> m_states;

    /// The chemical states in the chemical field
    std::vector<ChemicalProperties> m_properties;

    /// Construct a ChemicalField::Impl instance.
    Impl(Index size, const ChemicalSystem& system)
    : m_size(size), m_system(system), m_states(size, ChemicalState(system)), m_properties(size, ChemicalProperties(system))
    {
    }
    /// Construct a ChemicalField::Impl instance.
    Impl(Index size, const ChemicalState& state)
    : m_size(size), m_system(state.system()), m_states(size, state), m_properties(size, state.properties())
    {
    }

    auto set(const ChemicalState& state) -> void
    {
        for(auto& item : m_states)
            item = state;
    }

    auto temperature(VectorRef values) -> void
    {
        const Index len = m_size;
        for(Index i = 0; i < len; ++i)
            values[i] = m_states[i].temperature();
    }

    auto pressure(VectorRef values) -> void
    {
        const Index len = m_size;
        for(Index i = 0; i < len; ++i)
            values[i] = m_states[i].pressure();
    }

    auto elementAmounts(VectorRef values) -> void
    {
        const Index len = m_size;
        const Index num_elements = m_system.numElements();
        Index offset = 0;
        for(Index i = 0; i < len; ++i, offset += num_elements)
            values.segment(offset, num_elements) = m_states[i].elementAmounts();
    }

    auto output(const std::string& filename, const StringList& quantities) -> void
    {
        ChemicalOutput out(m_system);
        out.filename(filename);
        for(const auto& quantity : quantities)
            out.add(quantity);
    }

};

/// Construct a ChemicalField instance using Impl
ChemicalField::ChemicalField(Index size, const ChemicalSystem& system)
: pimpl(new Impl(size, system))
{
}

/// Construct a ChemicalField instance using Impl
ChemicalField::ChemicalField(Index size, const ChemicalState& state)
: pimpl(new Impl(size, state))
{
}
// Copy constructor
ChemicalField::ChemicalField(const ChemicalField& other)
: pimpl(new Impl(*other.pimpl))
{
}

ChemicalField::~ChemicalField()
{}

auto ChemicalField::operator=(ChemicalField other) -> ChemicalField&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalField::begin() const -> ConstIterator
{
    return pimpl->m_states.cbegin();
}

auto ChemicalField::begin() -> Iterator
{
    return pimpl->m_states.begin();
}

auto ChemicalField::end() const -> ConstIterator
{
    return  pimpl->m_states.cend();
}

auto ChemicalField::end() -> Iterator
{
    return pimpl->m_states.end();
}

auto ChemicalField::operator[](Index index) const -> const ChemicalState&
{
    return  pimpl->m_states[index];
}

auto ChemicalField::operator[](Index index) -> ChemicalState&
{
    return pimpl->m_states[index];
}

auto ChemicalField::size() const -> Index
{
    return pimpl->m_size;
}
auto ChemicalField::set(const ChemicalState& state) -> void
{
    pimpl->set(state);
}

auto ChemicalField::temperature(VectorRef values) -> void
{
    pimpl->temperature(values);
}

auto ChemicalField::pressure(VectorRef values) -> void
{
    pimpl->pressure(values);
}

auto ChemicalField::elementAmounts(VectorRef values) -> void
{
    pimpl->elementAmounts(values);
}

auto ChemicalField::output(const std::string& filename, const StringList& quantities) -> void
{
    pimpl->output(filename, quantities);
}

auto ChemicalField::properties() -> std::vector<ChemicalProperties>&
{
    return pimpl->m_properties;
}
auto ChemicalField::states() -> std::vector<ChemicalState>&
{
    return pimpl->m_states;
}

} // namespace Reaktoro
