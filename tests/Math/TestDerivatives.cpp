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

#include "TestDerivatives.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

#define ARMA_ASSERT_EQUAL_DELTA(expected, actual, delta) ASSERT(arma::norm(expected-actual)/arma::norm(expected) < delta)

} // namespace

auto test_derivativeForwardScalar() -> void
{
    ScalarFunction f = [](const Vector& x) { return arma::sum(arma::log(x)); };
    const Vector x = {1e-8, 1e8};
    const Vector expected = 1/x;
    const Vector actual = derivativeForward(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-6);
}

auto test_derivativeBackwardScalar() -> void
{
    ScalarFunction f = [](const Vector& x) { return arma::sum(arma::log(x)); };
    const Vector x = {1e-8, 1e8};
    const Vector expected = 1/x;
    const Vector actual = derivativeBackward(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-6);
}

auto test_derivativeCentralScalar() -> void
{
    ScalarFunction f = [](const Vector& x) { return arma::sum(arma::log(x)); };
    const Vector x = {1e-8, 1e8};
    const Vector expected = 1/x;
    const Vector actual = derivativeCentral(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-8);
}

auto test_derivativeForwardVector() -> void
{
    VectorFunction f = [](const Vector& x) { return arma::log(x); };
    const Vector x = {1e-8, 1e8};
    const Matrix expected = arma::diagmat(1/x);
    const Matrix actual = derivativeForward(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-6);
}

auto test_derivativeBackwardVector() -> void
{
    VectorFunction f = [](const Vector& x) { return arma::log(x); };
    const Vector x = {1e-8, 1e8};
    const Matrix expected = arma::diagmat(1/x);
    const Matrix actual = derivativeBackward(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-6);
}

auto test_derivativeCentralVector() -> void
{
    VectorFunction f = [](const Vector& x) { return arma::log(x); };
    const Vector x = {1e-8, 1e8};
    const Matrix expected = arma::diagmat(1/x);
    const Matrix actual = derivativeCentral(f, x);
    ARMA_ASSERT_EQUAL_DELTA(expected, actual, 1e-8);
}

auto testSuiteDerivatives() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_derivativeForwardScalar);
    s += CUTE(test_derivativeBackwardScalar);
    s += CUTE(test_derivativeCentralScalar);
    s += CUTE(test_derivativeForwardVector);
    s += CUTE(test_derivativeBackwardVector);
    s += CUTE(test_derivativeCentralVector);

    return s;
}

} // namespace Reaktoro
