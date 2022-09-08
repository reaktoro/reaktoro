// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "KineticsSensitivity.hpp"

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>

namespace Reaktoro {

KineticsSensitivity::KineticsSensitivity()
: EquilibriumSensitivity()
{}

KineticsSensitivity::KineticsSensitivity(EquilibriumSpecs const& specs)
: EquilibriumSensitivity(specs)
{}

KineticsSensitivity::KineticsSensitivity(EquilibriumSensitivity const& other)
: EquilibriumSensitivity(other)
{}

auto KineticsSensitivity::initialize(EquilibriumSpecs const& specs) -> void
{
    EquilibriumSensitivity::initialize(specs);
}

auto KineticsSensitivity::dnds(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::dnds() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticsSensitivity::dnddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::dpds(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::dpds() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticsSensitivity::dpddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::dqds(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::dqds() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticsSensitivity::dqddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticsSensitivity::duds() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticsSensitivity::duddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

} // namespace Reaktoro
