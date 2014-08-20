/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "PseudoInverse.hpp"

// Eigen includes
#include <Eigen/Dense>
using namespace Eigen;

namespace Reaktor {

auto pseudoInverse(const Matrix& M) -> Matrix
{
    JacobiSVD<Matrix> jacobi(M, ComputeThinU | ComputeThinV);

    const double epsilon = 1.0e-16;
    const double tol = epsilon * std::max(M.rows(), M.cols()) * jacobi.singularValues().cwiseAbs().maxCoeff();

    return jacobi.matrixV() * (jacobi.singularValues().array() > tol).
        select(jacobi.singularValues().cwiseInverse(), 0.0).asDiagonal() * jacobi.matrixU().transpose();
}

} /* namespace Optima */
