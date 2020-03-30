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
#include <deque>
#include <optional>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>

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
    /// Construct a default SubstanceCriticalProps instance.
    SubstanceCriticalProps();

    /// Construct a SubstanceCriticalProps instance with given name and chemical formula.
    SubstanceCriticalProps(std::string name, const ChemicalFormula& formula);

    /// Construct a SubstanceCriticalProps instance with given name and chemical formula.
    SubstanceCriticalProps(std::string name, const ChemicalFormula& formula, const SubstanceCriticalPropsData& data);

    /// Set the critical temperature of the substance (in K)
    auto setTemperature(real value) -> void;

    /// Set the critical temperature of the substance with given unit.
    auto setTemperature(real value, std::string unit) -> void;

    /// Set the critical pressure of the substance (in Pa)
    auto setPressure(real value) -> void;

    /// Set the critical pressure of the substance with given unit.
    auto setPressure(real value, std::string unit) -> void;

    /// Set the acentric factor of the substance.
    auto setAcentricFactor(real value) -> void;

    /// Return the name of the substance.
    auto name() const -> const std::string&;

    /// Return the chemical formula of the substance.
    auto formula() const -> const ChemicalFormula&;

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
    /// The name of the substance.
    std::string m_name;

    /// The chemical formula of the substance.
    ChemicalFormula m_formula;

    /// The critical properties data of the substance.
    SubstanceCriticalPropsData m_data;
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

    /// Return the number of substances in the database.
    static auto size() -> std::size_t;

    /// Return the substance and its critical properties with given name.
    static auto getWithName(const std::string& name) -> std::optional<SubstanceCriticalProps>;

    /// Return the substance and its critical properties with given chemical formula.
    static auto getWithFormula(const ChemicalFormula& formula) -> std::optional<SubstanceCriticalProps>;

    /// Return the substance and its critical properties with given name or chemical formula.
    static auto get(const std::string& name_or_formula) -> std::optional<SubstanceCriticalProps>;

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
