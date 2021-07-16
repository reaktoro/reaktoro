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
#include <deque>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// A type used to represent the critical properties of a substance.
struct SubstanceCriticalPropsData
{
    /// The critical temperature of the substance (in K).
    real Tcr = {};

    /// The critical pressure of the substance (in Pa).
    real Pcr = {};

    /// The acentric factor of the substance.
    real omega = {};
};

/// A type used to represent a substance and its critical properties.
class SubstanceCriticalProps
{
public:
    /// Construct a SubstanceCriticalProps instance.
    SubstanceCriticalProps(const StringList& names);

    /// Construct a SubstanceCriticalProps instance with given data.
    /// The given names will be converted to uppercase, suffix will be removed,
    /// and spaces will be replaced by dashes. So, for example, the substance name
    /// `carbon dioxide` is replaced by `CARBON-DIOXIDE` and `HCl(g)` by `HCL`.
    /// @param data The critical property data of the substance
    /// @param names The names that can uniquely identify the substance *(case-insensitive)*
    SubstanceCriticalProps(const SubstanceCriticalPropsData& data, const StringList& names);

    /// Set the critical temperature of the substance (in K)
    auto setTemperature(real value) -> void;

    /// Set the critical temperature of the substance with given unit.
    auto setTemperature(real value, String unit) -> void;

    /// Set the critical pressure of the substance (in Pa)
    auto setPressure(real value) -> void;

    /// Set the critical pressure of the substance with given unit.
    auto setPressure(real value, String unit) -> void;

    /// Set the acentric factor of the substance.
    auto setAcentricFactor(real value) -> void;

    /// Return the names that uniquely identify the substance.
    auto names() const -> const Strings&;

    /// Return the critical temperature of the substance (in K)
    auto temperature() const -> real;

    /// Return the critical pressure of the substance (in Pa)
    auto pressure() const -> real;

    /// Return the acentric factor of the substance.
    auto acentricFactor() const -> real;

    /// Return the critical properties data of the substance.
    auto data() const -> const SubstanceCriticalPropsData&;

    /// Return the critical properties data of the substance.
    operator SubstanceCriticalPropsData() const;

private:
    /// The critical properties data of the substance.
    SubstanceCriticalPropsData m_data;

    /// The names that uniquely identify the substance.
    Strings m_names;
};

/// A type used store a collection of substances and their critical properties.
class CriticalProps
{
public:
    /// Construct a copy of a CriticalProps object [deleted].
    CriticalProps(const CriticalProps&) = delete;

    /// Assign a CriticalProps object to this [deleted].
    auto operator=(const CriticalProps&) -> CriticalProps& = delete;

    /// Return the single CriticalProps object.
    static auto instance() -> CriticalProps&;

    /// Return the substances and their critical properties in the database.
    static auto substances() -> const std::deque<SubstanceCriticalProps>&;

    /// Append a substance and its critical properties in to the database.
    static auto append(SubstanceCriticalProps substance) -> void;

    /// Append a substance and its critical properties in the database or overwrite an existing one with conflicting name.
    static auto overwrite(SubstanceCriticalProps substance) -> void;

    /// Return the number of substances in the database.
    static auto size() -> Index;

    /// Return the substance and its critical properties with given name (e.g. "WATER", "CARBON-DIOXIDE", "HYDROGEN-SULFIDE", etc.).
    static auto get(String name) -> Optional<SubstanceCriticalProps>;

    /// Return the substance and its critical properties with given alternative names.
    static auto get(const StringList& names) -> Optional<SubstanceCriticalProps>;

    /// Return begin const iterator of this CriticalProps instance
    auto begin() const;

    /// Return begin iterator of this CriticalProps instance
    auto begin();

    /// Return end const iterator of this CriticalProps instance
    auto end() const;

    /// Return end iterator of this CriticalProps instance
    auto end();

private:
    /// The substances stored in the database.
    std::deque<SubstanceCriticalProps> m_substances;

private:
    /// Construct a default CriticalProps object [private].
    CriticalProps();

    /// Destroy this CriticalProps object [private].
    ~CriticalProps();
};

} // namespace Reaktoro
