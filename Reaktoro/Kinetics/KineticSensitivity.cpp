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

#include "KineticSensitivity.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

struct KineticSensitivity::Impl
{
    /// Construct a default Impl object.
    Impl()
    {}

    /// Construct a default Impl object.
    Impl(EquilibriumSpecs const& specs)
    {
        initialize(specs);
    }

    /// Initialize this KineticSensitivity object with given equilibrium problem specifications.
    auto initialize(EquilibriumSpecs const& specs) -> void
    {

    }
};

KineticSensitivity::KineticSensitivity()
: pimpl(new Impl())
{}

KineticSensitivity::KineticSensitivity(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

KineticSensitivity::KineticSensitivity(KineticSensitivity const& other)
: pimpl(new Impl(*other.pimpl))
{}

KineticSensitivity::~KineticSensitivity()
{}

auto KineticSensitivity::operator=(KineticSensitivity other) -> KineticSensitivity&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticSensitivity::initialize(EquilibriumSpecs const& specs) -> void
{
    pimpl->initialize(specs);
}

auto KineticSensitivity::dndS(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dndS() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticSensitivity::dnddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dpdS(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dpdS() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticSensitivity::dpddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dqdS(String const& Sid) const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dqdS() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticSensitivity::dqddt() const -> VectorXdConstRef
{
    return VectorXd{};
}

auto KineticSensitivity::dudS() const -> MatrixXdConstRef
{
    return MatrixXd{};
}

auto KineticSensitivity::duddt() const -> VectorXdConstRef
{
    return VectorXd{};
}


} // namespace Reaktoro
