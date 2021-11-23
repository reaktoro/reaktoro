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

/// A type used to describe the state of an complexation surface.
/// @see ComplexationComposition
struct ComplexationSurfaceState
{
    /// The equivalent fractions of the exchange species
    ArrayXr beta;

    /// The natural logarithms of the activity coefficients (calculated during the activity model evaluation)
    ArrayXr lng;
};

/// A type used to describe an complexation surface.
/// The ComplexationSurface class is defined as a collection of Species objects, representing,
/// therefore, a composition of complexation phase. Its main purpose is to provide the
/// necessary operations in the calculation of activities of complexation species.
/// It implements methods for the calculation of equivalent fractions of species in ionic
/// exchange phase. In addition, it provides methods that retrieves information about the
/// exchanger (e.g., X-) and exchange species (e.g., NaX, CaX2).
class ComplexationSurface
{
public:
    /// Construct a default ComplexationSurface instance.
    ComplexationSurface();

    /// Construct an ComplexationSurface instance with given species.
    explicit ComplexationSurface(const SpeciesList& species);

    /// Return a deep copy of this ComplexationSurface object.
    auto clone() const -> ComplexationSurface;

    /// Return the exchange species on the surface with given index.
    /// @param idx The index of the species in the complexation surface
    auto species(Index idx) const -> const Species&;

    /// Return the exchange species on the surface.
    auto species() const -> const SpeciesList&;

    /// Return the array of exchanger's equivalents numbers (or cation charges) in complexation species.
    auto ze() const -> ArrayXdConstRef;

    /// Calculate the state of the aqueous mixture.
    /// @param T The temperature (in K)
    /// @param P The pressure (in Pa)
    /// @param x The fraction of the species in the composition
    auto state(real T, real P, ArrayXrConstRef x) -> ComplexationSurfaceState;

    /// Return the state of the aqueous mixture.
    auto state() const -> ComplexationSurfaceState;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro
