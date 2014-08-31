/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class EquilibriumConstraints;
class Partitioning;

class EquilibriumProblem
{
public:
    explicit EquilibriumProblem(const ChemicalSystem& system);

    EquilibriumProblem(const ChemicalSystem& system, const Partitioning& partitioning);

    EquilibriumProblem(const EquilibriumProblem& other);

    ~EquilibriumProblem();

    auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

    auto addCompound(std::string compound, units::Amount value) -> void;

    auto addCompound(std::string compound, units::Mass value) -> void;

    auto addCompound(std::string compound, double value, std::string unit) -> void;

    auto addSpecies(Index idx, units::Amount value) -> void;

    auto addSpecies(Index idx, units::Mass value) -> void;

    auto addSpecies(std::string species, units::Amount value) -> void;

    auto addSpecies(std::string species, units::Mass value) -> void;

    auto addSpecies(std::string species, double value, std::string unit) -> void;

    auto add(const ChemicalState& state) -> void;

    auto add(std::string entity, double value, std::string unit) -> void;

    auto setAcidity(double value) -> void;

    auto setChargeBalance() -> void;

    auto setChargeBalance(std::string phase) -> void;

    auto setSpecies(std::string species, units::Mass value) -> void;

    auto setSpecies(std::string species, units::Amount value) -> void;

    auto setActivity(std::string species, double value) -> void;

    auto setPartialPressure(std::string gas, units::Pressure value) -> void;

    auto setFreeElements(std::string elements) -> void;

    auto setFreeElements(std::vector<std::string> elements) -> void;

    auto freeElements() const -> const std::vector<std::string>&;

    auto infoConstraints() -> std::vector<std::string>;

    auto constraints() const -> EquilibriumConstraints;

    operator EquilibriumConstraints() const;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
