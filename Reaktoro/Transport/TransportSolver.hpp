// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;

class BoundaryState
{
public:
    BoundaryState(const ChemicalState& state);

private:
};

class BoundaryFlux
{
public:
    BoundaryFlux(const ChemicalState& state);

private:
};

class ChemicalField
{
public:
    ChemicalField(Index size);

    ChemicalField(Index size, const ChemicalState& state);

    auto temperatures(VectorRef values) -> void;

    auto pressures(VectorRef values) -> void;

    auto elementAmounts(MatrixRef values) -> void;

    auto states() const -> std::vector<std::unique_ptr<const ChemicalState>>;

    auto states() -> std::vector<std::unique_ptr<ChemicalState>>;

    auto state(Index index) const -> const ChemicalState&;

    auto state(Index index) -> ChemicalState&;

private:
};

} // namespace Reaktoro
