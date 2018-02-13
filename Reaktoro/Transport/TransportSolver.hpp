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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
//class ChemicalProperties;
//class ChemicalState;
//class ChemicalSystem;

//class BoundaryState
//{
//public:
//    BoundaryState(const ChemicalState& state);
//
//private:
//};
//
//class BoundaryFlux
//{
//public:
//    BoundaryFlux(const ChemicalState& state);
//
//private:
//};


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

private:
    /// The number of degrees of freedom in the chemical field.
    Index m_size;

    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem m_system;

    /// The chemical states in the chemical field
    std::vector<ChemicalState> m_states;

    /// The chemical states in the chemical field
    std::vector<ChemicalProperties> m_properties;
};

class TridiagonalMatrix
{
public:
    TridiagonalMatrix(Index size) : m_size(size), m_data(size * 3) {}

    auto size() const -> Index { return m_size; }

    auto data() -> VectorRef { return m_data; }

    auto data() const -> VectorConstRef { return m_data; }

    auto row(Index index) -> VectorRef { return m_data.segment(3 * index, 3); }

    auto row(Index index) const -> VectorConstRef { return m_data.segment(3 * index, 3); }

    auto a() -> VectorStridedRef { return Vector::Map(m_data.data() + 3, m_size - 1, Eigen::InnerStride<3>()); }

    auto a() const -> VectorConstRef { return Vector::Map(m_data.data() + 3, m_size - 1, Eigen::InnerStride<3>()); }

    auto b() -> VectorStridedRef { return Vector::Map(m_data.data() + 1, m_size, Eigen::InnerStride<3>()); }

    auto b() const -> VectorConstRef { return Vector::Map(m_data.data() + 1, m_size, Eigen::InnerStride<3>()); }

    auto c() -> VectorStridedRef { return Vector::Map(m_data.data() + 2, m_size - 1, Eigen::InnerStride<3>()); }

    auto c() const -> VectorConstRef { return Vector::Map(m_data.data() + 2, m_size - 1, Eigen::InnerStride<3>()); }

    auto factorize() -> void;

    auto solve(VectorRef x, VectorConstRef b) const -> void;

    auto solve(VectorRef x) const -> void;

    operator Matrix() const;

private:
    Index m_size;

    Vector m_data;
};

///
class TransportSolver
{
public:
    /// Construct a default TransportSolver instance.
    TransportSolver(const ChemicalSystem& system);

    auto setVelocity(double val) -> void { m_velocity = val; }

    auto setDiffusionCoeff(double val) -> void { m_diffusion = val; }

    auto setBoundaryState(const ChemicalState& state) -> void;

    auto setTimeStep(double val) -> void { m_dt = val; }

    auto initialize(const ChemicalField& field) -> void;

    auto step(ChemicalField& field) -> void;

private:
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem m_system;

    /// The velocity in the transport problem (in m/s).
    double m_velocity = 0.0;

    /// The diffusion coefficient in the transport problem (in m^2/s).
    double m_diffusion = 0.0;

    /// The time step used to solve the transport problem (in s).
    double m_dt = 0.0;

    /// The chemical state on the boundary.
    ChemicalState m_bc;

    /// The amounts of elements on the boundary.
    Vector m_bbc;

    /// The amounts of elements on each cell of the mesh (each column corresponds to a mesh cell).
    Matrix m_b;
};

} // namespace Reaktoro
