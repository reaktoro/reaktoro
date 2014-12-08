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
namespace {

Indices linearlyIndependentCols(const Eigen::MatrixXd& A)
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);
    Indices indices(I.data(), I.data() + rank);
    return indices;
}

} // namespace

Indices linearlyIndependentCols(const Matrix& A)
{
    Eigen::MatrixXd A_eigen = Eigen::Map<const Eigen::MatrixXd>(A.memptr(), A.n_rows, A.n_cols);
    return linearlyIndependentCols(A_eigen);
}

Indices linearlyIndependentRows(const Matrix& A)
{
    Eigen::MatrixXd At_eigen = Eigen::Map<const Eigen::MatrixXd>(A.memptr(), A.n_rows, A.n_cols).transpose();
    return linearlyIndependentCols(At_eigen);
}

Indices linearlyIndependentCols(const Matrix& A, Matrix& B)
{
    Indices indices = linearlyIndependentCols(A);
    B.resize(A.n_rows, indices.size());
    for(unsigned i = 0; i < B.n_cols; ++i)
        B.col(i) = A.col(indices[i]);
    return indices;
}

Indices linearlyIndependentRows(const Matrix& A, Matrix& B)
{
    Indices indices = linearlyIndependentRows(A);
    B.resize(indices.size(), A.n_cols);
    for(unsigned i = 0; i < B.n_rows; ++i)
        B.row(i) = A.row(indices[i]);
    return indices;
}

} // namespace Reaktor
