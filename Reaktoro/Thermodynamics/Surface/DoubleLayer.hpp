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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used to describe the state of an ion exchange surface.
/// @see IonExchangeComposition
struct DoubleLayerState
{
    /// The temperature of the DDL (in K).
    real T;

    /// The pressure of the DDL (in Pa).
    real P;

    /// The effective ionic strength of the aqueous mixture (in mol/kg)
    real Ie;

    /// The stoichiometric ionic strength of the aqueous mixture (in mol/kg)
    real Is;

    /// The molalities of the aqueous species in DDL (in mol/kg)
    ArrayXr m;

    /// The stoichiometric molalities of the ionic species in DDL (in mol/kg)
    ArrayXr ms;

    /// The activities of the species in double layer
    ArrayXr a;

    /// The natural logarithms of the activity coefficients
    ArrayXr lng;
};

/// A type used to describe an diffusive double layer between complexation surface and aqueous solution.
/// The DoubleLayer class is defined as a collection of Species objects, representing,
/// therefore, a composition of the DoubleLayer phase. Its main purpose is to provide the
/// necessary operations in the calculation of activities of diffusive double layer species.
class DoubleLayer
{
public:
    /// Construct a default DoubleLayer instance.
    DoubleLayer();

    /// Construct an DoubleLayer instance with given species.
    explicit DoubleLayer(const SpeciesList& species);

    /// Return a deep copy of this DoubleLayer object.
    auto clone() const -> DoubleLayer;

    /// Return the exchange species on the surface with given index.
    /// @param idx The index of the species in the DDL
    auto species(Index idx) const -> const Species&;

    /// Return the exchange species on the surface.
    auto species() const -> const SpeciesList&;

    /// Return the array of DDL species' charges.
    auto charges() const -> ArrayXdConstRef;

    /// Return the index of the solvent species in the mixture.
    auto indexWater() const -> Index;

    /// Calculate the state of the aqueous mixture.
    /// @param T The temperature (in K)
    /// @param P The pressure (in Pa)
    /// @param x The fraction of the species in the composition
    auto state(real T, real P, ArrayXrConstRef x) -> DoubleLayerState;

    /// Return the state of the aqueous mixture.
    auto state() const -> DoubleLayerState;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro
