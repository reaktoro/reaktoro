// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "MathUtils.hpp"

// Eigen includes
#include <Reaktor/eigen/Dense>

namespace Reaktor {

auto linearlyIndependentCols(const Matrix& A) -> Indices
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);
    Indices indices(I.data(), I.data() + rank);
    return indices;
}

auto linearlyIndependentRows(const Matrix& A) -> Indices
{
    const Matrix At = A.transpose();
    return linearlyIndependentCols(At);
}

auto linearlyIndependentCols(const Matrix& A, Matrix& B) -> Indices
{
    Indices indices = linearlyIndependentCols(A);
    Matrix C(A.rows(), indices.size());
    for(unsigned i = 0; i < indices.size(); ++i)
        C.col(i) = A.col(indices[i]);
    B.noalias() = C;
    return indices;
}

auto linearlyIndependentRows(const Matrix& A, Matrix& B) -> Indices
{
    Indices indices = linearlyIndependentRows(A);
    Matrix C(indices.size(), A.cols());
    for(unsigned i = 0; i < indices.size(); ++i)
        C.row(i) = A.row(indices[i]);
    B.noalias() = C;
    return indices;
}

auto inverseShermanMorrison(const Matrix& invA, const Vector& D) -> Matrix
{
    Matrix invM = invA;
    for(unsigned i = 0; i < D.rows(); ++i)
        invM = invM - (D[i]/(1 + D[i]*invM(i, i)))*invM.col(i)*invM.row(i);
    return invM;
}

} // namespace Reaktor
