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

// C++ includes
#include <functional>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// Compute the largest step length @f$\alpha@f$ such that
/// @f$\mathbf{p} + \alpha\Delta\mathbf{p}@f$ is on the
/// lower bound @f$\mathbf{x}_l=\mathbf{0}@f$.
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
auto largestStep(VectorXdConstRef p, VectorXdConstRef dp) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\mathbf{p}+\alpha\Delta\mathbf{p}\geq(1-\tau)\mathbf{p}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, double tau) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\mathbf{p}+\alpha\Delta\mathbf{p}\geq(1-\tau)\mathbf{p}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
/// @param ilimiting The index of the limiting variable
auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, double tau, Index& ilimiting) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\alpha C\Delta\mathbf{p}\geq-\tau C\mathbf{p}+\mathbf{r}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param C The left-hand side matrix that defines the inequality constraint @f$C\mathbf{p}\geq\mathbf{r}@f$
/// @param r The right-hand side vector that defines the inequality constraint @f$C\mathbf{p}\geq\mathbf{r}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
auto fractionToTheBoundary(VectorXdConstRef p, VectorXdConstRef dp, MatrixXdConstRef C, VectorXdConstRef r, double tau) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\mathbf{p}+\alpha\Delta\mathbf{p}\geq(1-\tau)\mathbf{p}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param lower The lower bound for the step @f$\Delta\mathbf{p}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
auto fractionToTheLowerBoundary(VectorXdConstRef p, VectorXdConstRef dp, VectorXdConstRef lower, double tau) -> double;

/// Check if a float number is less than another by a base value.
/// This method is particularly useful when comparing two float numbers
/// where round-off errors can compromise the checking.
/// The following is used for the comparison:
/// @f[a < b + 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto lessThan(double a, double b, double baseval) -> bool;

/// Check if a float number is greater than another by a base value.
/// This method is particularly useful when comparing two float numbers
/// where round-off errors can compromise the checking.
/// The following is used for the comparison:
/// @f[a > b - 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto greaterThan(double a, double b, double baseval) -> bool;

/// Return the floating-point representation of positive infinity
auto infinity() -> double;

/// Return an inverse Hessian function based on the BFGS Hessian approximation
auto bfgs() -> std::function<MatrixXd(VectorXdConstRef, VectorXdConstRef)>;

/// Calculate the minimum of a single variable function using the Golden Section Search algorithm.
auto minimizeGoldenSectionSearch(const std::function<double(double)>& f, double a, double b, double tol = 1e-5) -> double;

/// Calculate the minimum of a single variable function using the Brent algorithm.
auto minimizeBrent(const std::function<double(double)>& f, double min, double max, double tolerance = 1e-5, unsigned maxiters = 100) -> double;

} // namespace Reaktoro
