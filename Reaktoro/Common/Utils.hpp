// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

namespace Reaktoro {

/// Return the mole fractions of the species with given amounts.
inline auto moleFractions(ArrayXrConstRef n) -> ArrayXr
{
    const auto nspecies = n.size();
    if(nspecies == 1)
        return ArrayXr{{1.0}};
    const auto nsum = n.sum();
    if(nsum != 0.0) return n/nsum;
    else return ArrayXr::Zero(nspecies);
}

/// Return the derivatives of the mole fractions of the species with respect to their amounts.
inline auto moleFractionsJacobian(ArrayXrRef n) -> MatrixXd
{
    using autodiff::jacobian;
    using autodiff::wrt;
    using autodiff::at;
    return jacobian(moleFractions, wrt(n), at(n));
}

} // namespace Reaktoro
