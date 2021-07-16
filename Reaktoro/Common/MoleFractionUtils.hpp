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
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// Compute the mole fractions of the species with given species amounts.
/// @param n The vector with the species amounts.
/// @param[out] x The vector of species mole fractions.
template<typename ArrayConstRef, typename ArrayRef>
auto moleFractions(ArrayConstRef&& n, ArrayRef&& x)
{
    const auto nspecies = n.size();
    if(nspecies == 1) {
        x.fill(1.0);
        return;
    }
    const auto nsum = n.sum();
    if(nsum != 0.0) x = n/nsum;
    else x.fill(0.0);
}

/// Compute the mole fractions of the species with given species amounts.
/// @param n The vector with the species amounts.
template<typename ArrayConstRef>
auto moleFractions(ArrayConstRef&& n)
{
    using T = Decay<decltype(n[0])>;
    const auto N = n.size();
    ArrayX<T> x(N);
    moleFractions(n, x);
    return x;
}

/// Compute the Jacobian matrix of the species mole fractions (@eq{J=\frac{\partial x}{\partial n}}).
/// @param n The vector with the species amounts.
/// @param[out] J The output Jacobian matrix.
template<typename ArrayConstRef, typename MatrixRef>
auto moleFractionsJacobian(ArrayConstRef&& n, MatrixRef&& J)
{
    //-----------------------------------------------------------------------------------------------------
    // == LATEX ==
    //-----------------------------------------------------------------------------------------------------
    // \frac{\partial x_{i}}{\partial n_{j}}=\frac{x_{i}}{n_{i}}\left(\delta_{ij}-x_{i}\right)=\frac{1}{n_{\Sigma}}\left(\delta_{ij}-x_{i}\right)
    //-----------------------------------------------------------------------------------------------------
    const auto N = n.size();
    assert(J.rows() == N);
    assert(J.cols() == N);
    const auto nsum = n.sum();
    if(nsum == 0.0) J.fill(0.0);
    else for(auto i = 0; i < N; ++i) {
        const auto ni = n[i];
        const auto xi = ni/nsum;
        const auto aux = -xi/nsum;
        for(auto j = 0; j < N; ++j)
            J(i, j) = aux; // J(i, j) = -xi/nsum
        J(i, i) = (1 - xi)/nsum; // J(i, i) = (1 - xi)/nsum
    }
}

/// Compute the Jacobian matrix of the species mole fractions (@eq{J=\frac{\partial x}{\partial n}}).
/// @param n The vector with the species amounts.
template<typename ArrayConstRef>
auto moleFractionsJacobian(ArrayConstRef&& n)
{
    using T = Decay<decltype(n[0])>;
    const auto N = n.size();
    MatrixX<T> J(N, N);
    moleFractionsJacobian(n, J);
    return J;
}

/// Compute the Jacobian matrix of the species mole fractions in natural log (@eq{J=\frac{\partial\ln x}{\partial n}}).
/// @param n The vector with the species amounts.
/// @param[out] J The output Jacobian matrix.
template<typename ArrayConstRef, typename MatrixRef>
auto lnMoleFractionsJacobian(ArrayConstRef&& n, MatrixRef&& J) -> void
{
    //-----------------------------------------------------------------------------------------------------
    // == LATEX ==
    //-----------------------------------------------------------------------------------------------------
    // \frac{\partial\ln x_{i}}{\partial n_{j}}=\frac{\delta_{ij}}{n_{i}}-\frac{1}{n_{\Sigma}}=\frac{1}{n_{i}}\left(\delta_{ij}-x_{i}\right)
    //-----------------------------------------------------------------------------------------------------
    using T = Decay<decltype(J(0, 0))>;
    const auto N = n.size();
    assert(J.rows() == N);
    assert(J.cols() == N);
    const auto nsum = n.sum();
    if(nsum == 0.0) J.fill(0.0);
    else for(auto i = 0; i < N; ++i) {
        const auto ni = n[i];
        const auto xi = ni/nsum;
        const auto aux = -xi/ni;
        for(auto j = 0; j < N; ++j)
            J(i, j) = static_cast<T>(aux); // J(i, j) = -xi/ni
        J(i, i) = static_cast<T>((1 - xi)/ni); // J(i, i) = (1 - xi)/ni
    }
}

/// Compute the Jacobian matrix of the species mole fractions in natural log (@eq{J=\frac{\partial\ln x}{\partial n}}).
/// @param n The vector with the species amounts.
template<typename ArrayConstRef>
auto lnMoleFractionsJacobian(ArrayConstRef&& n)
{
    using T = Decay<decltype(n[0])>;
    const auto N = n.size();
    MatrixX<T> J(N, N);
    lnMoleFractionsJacobian(n, J);
    return J;
}

} // namespace Reaktoro
