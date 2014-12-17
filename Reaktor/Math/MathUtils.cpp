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

Indices linearlyIndependentCols(const Matrix& A)
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);
    Indices indices(I.data(), I.data() + rank);
    return indices;
}

Indices linearlyIndependentRows(const Matrix& A)
{
    const Matrix At = A.transpose();
    return linearlyIndependentCols(At);
}

Indices linearlyIndependentCols(const Matrix& A, Matrix& B)
{
    Indices indices = linearlyIndependentCols(A);
    B.resize(A.rows(), indices.size());
    for(unsigned i = 0; i < B.cols(); ++i)
        B.col(i) = A.col(indices[i]);
    return indices;
}

Indices linearlyIndependentRows(const Matrix& A, Matrix& B)
{
    Indices indices = linearlyIndependentRows(A);
    B.resize(indices.size(), A.cols());
    for(unsigned i = 0; i < B.rows(); ++i)
        B.row(i) = A.row(indices[i]);
    return indices;
}

} // namespace Reaktor
