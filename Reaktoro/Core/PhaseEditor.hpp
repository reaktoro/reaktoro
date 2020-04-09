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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/ThermoEngine.hpp>

namespace Reaktoro {

/// The auxiliary type used to specify phase species to be determined from element symbols.
struct Speciate
{
    /// The symbols of the elements composing the species in a phase.
    Strings symbols;
};

/// The auxiliary function used to specify phase species to be determined from element symbols.
inline auto speciate(const StringList& symbols) { return Speciate{symbols}; };

/// The base type for all other classes defining phase types.
/// @ingroup Core
class PhaseEditor
{
public:
    /// Construct a PhaseEditor object.
    explicit PhaseEditor(const StringList& species);

    /// Construct a PhaseEditor object.
    explicit PhaseEditor(const Speciate& elements);

    /// Destroy this PhaseEditor object.
    virtual ~PhaseEditor();

    /// Set a unique name for the phase.
    auto name(String name) -> PhaseEditor&;

    /// Set the state of matter for the phase.
    auto set(StateOfMatter state) -> PhaseEditor&;

    /// Set the activity model for the phase.
    auto set(const ActivityModel& model) -> PhaseEditor&;

    /// Return the name of the phase.
    auto name() const -> String;

    /// Return the state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the names of the selected species to compose the phase (empty if not given).
    auto species() const -> const Strings&;

    /// Return the element symbols for automatic species selection (empty if not given).
    auto elements() const -> const Strings&;

    /// Return the specified activity model of the phase.
    auto activityModel() const -> ActivityModelPtr;

    /// Convert this PhaseEditor object into a Phase object.
    virtual auto convert(const ThermoEngine& engine, const Strings& elements) -> Phase;

private:
    /// The name of the phase.
    String m_name;

    /// The state of matter of the phase.
    StateOfMatter m_state = StateOfMatter::Solid;

    /// The names of the selected species to compose the phase.
    Strings m_species;

    /// The element symbols for automatic selection of the species composing the phase.
    Strings m_elements;

    /// The activity model of the phase.
    ActivityModelPtr m_activity_model_ptr;
};

} // namespace Reaktoro
