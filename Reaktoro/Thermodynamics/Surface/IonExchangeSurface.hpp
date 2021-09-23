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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used to describe an ion exchange surface.
/// The IonExchangeSurface class is defined as a collection of Species objects, representing,
/// therefore, a composition of ion exchange phase. Its main purpose is to provide the
/// necessary operations in the calculation of activities of ion exchange species.
/// It implements methods for the calculation of equivalence fractions of species in ionic
/// exchange phase. In addition, it provides methods that retrieves information about the
/// exchanger (e.g., X-) and exchange species (e.g., NaX, CaX2).
class IonExchangeSurface
{
public:
    /// Construct a default IonExchangeSurface instance.
    IonExchangeSurface();

    /// Construct an IonExchangeSurface instance with given species.
    explicit IonExchangeSurface(const SpeciesList& species);

    /// Return a deep copy of this IonExchangeSurface object.
    auto clone() const -> IonExchangeSurface;

    /// Return the exchange species on the surface with given index.
    auto species(Index idx) const -> const Species&;

    /// Return the exchange species on the surface.
    auto species() const -> const SpeciesList&;

    /// Return the indices of the exchanger on the surface.
    auto indexExchanger() const -> Index;

    /// Return the indices of the exchange species on the surface.
    auto indicesExchange() const -> const Indices&;

    /// Return the array of exchanger's equivalents numbers (or cation charges) in ion exchange species.
    auto ze() const -> ArrayXdConstRef;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro
