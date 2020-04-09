// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "PhaseEditor.hpp"

// Reaktoro includes

namespace Reaktoro {

PhaseEditor::PhaseEditor(const StringList& species)
: m_species(species)
{}

PhaseEditor::PhaseEditor(const Speciate& elements)
: m_elements(elements.symbols)
{}

PhaseEditor::~PhaseEditor()
{

}

auto PhaseEditor::name(String name) -> PhaseEditor&
{
    m_name = name;
    return *this;
}

auto PhaseEditor::set(StateOfMatter state) -> PhaseEditor&
{
    m_state = state;
    return *this;
}

auto PhaseEditor::set(const ActivityModel& model) -> PhaseEditor&
{
    m_activity_model_ptr = std::make_uniqmodel;
    return *this;

}

auto PhaseEditor::name() const -> String
{

}

auto PhaseEditor::stateOfMatter() const -> StateOfMatter
{

}

auto PhaseEditor::species() const -> const Strings&
{

}

auto PhaseEditor::elements() const -> const Strings&
{

}

auto PhaseEditor::activityPropsFn() const -> ActivityModelPtr
{

}

virtual PhaseEditor::auto convert(const ThermoEngine& engine, const Strings& elements) -> Phase
{

}

} // namespace Reaktoro
