// Reaktoro is a C++ library for computational reaction modelling.
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
#include <Reaktoro/Eigen/Dense>

namespace Reaktoro {

auto linearlyIndependentCols(const Matrix& A) -> Indices
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);
    Indices indices(I.data(), I.data() + rank);
    return indices;
}

auto linearlyIndependentColsPartialPivLU(const Matrix& A) -> Indices
{
    Matrix L, U;
    PermutationMatrix P, Q;
    lu(A);
    const Index n = A.cols();
    const Index m = A.rows();
    Indices indices;
    for(Index i = 0; i < m; ++i)
        for(Index j = i; j < n; ++j)
            if(U(i, j) != 0.0)
                { indices.push_back(j); break; }
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

/// Compute LU decomposition of a square matrix.
template<typename MatrixType>
auto luSquare(const MatrixType& A) -> DecompositionLU
{
    const Index m = A.rows();
    Eigen::PartialPivLU<Matrix> lu = A.lu();

    DecompositionLU res;
    res.L = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitLower>();
    res.U = lu.matrixLU().triangularView<Eigen::Upper>();
    res.P = lu.permutationP();
    res.Q.setIdentity(m);
    res.rank = m; // todo a proper test needs to be done to determine rank
    return res;
}

/// Compute LU decomposition of a rectangular matrix (more rows than cols).
template<typename MatrixType>
auto luRectangular1(const MatrixType& A) -> DecompositionLU
{
    const Index m = A.rows();
    const Index n = A.cols();
    Matrix M(n, n);
    M.topRows(m) = A;
    M.bottomLeftCorner(n-m, m).fill(0.0);
    M.bottomRightCorner(n-m, n-m).setIdentity();
    Eigen::PartialPivLU<Matrix> lu = M.lu();

    DecompositionLU res;
    res.L = lu.matrixLU().topLeftCorner(m, m).triangularView<Eigen::UnitLower>();
    res.U = lu.matrixLU().topRows(m).triangularView<Eigen::Upper>();
    res.P.resize(m);
    std::copy(lu.permutationP().indices().data(),
              lu.permutationP().indices().data() + m,
              res.P.indices().data());
    res.Q.setIdentity(n);
    for(Index i = 0; i < m; ++i)
        for(Index j = i; j < n; ++j)
            if(res.U(i, j) != 0.0)
            {
                auto tmp = res.Q.indices()[i];
                res.Q.indices()[i] = res.Q.indices()[j];
                res.Q.indices()[j] = tmp;
                res.rank = 1+i;
                break;
            }
    res.U = res.U * res.Q;
    return res;
}

/// Compute LU decomposition of a rectangular matrix (more cols than rows).
template<typename MatrixType>
auto luRectangular2(const MatrixType& A) -> DecompositionLU
{
    const Index m = A.rows();
    const Index n = A.cols();
    Matrix M(m, m);
    M.leftCols(n) = A;
    M.topRightCorner(n, m-n).fill(0.0);
    M.bottomRightCorner(m-n, m-n).setIdentity();
    Eigen::PartialPivLU<Matrix> lu = M.lu();

    DecompositionLU res;
    res.L = lu.matrixLU().leftCols(n).triangularView<Eigen::UnitLower>();
    res.U = lu.matrixLU().topLeftCorner(n, n).triangularView<Eigen::Upper>();
    res.P = lu.permutationP();
    res.Q.setIdentity(n);
    res.rank = n; // todo a proper test needs to be done to determine rank
    return res;
}

auto lu(const Matrix& A) -> DecompositionLU
{
    const Index m = A.rows();
    const Index n = A.cols();
    if(m == n) return luSquare(A);
    else if(m < n) return luRectangular1(A);
    else return luRectangular2(A);
}

auto fraction(double f, long maxden, long& num, long& den) -> void
{
    // Adapted from http://rosettacode.org/wiki/Convert_decimal_number_to_rational#C
    /*  a: continued fraction coefficients. */
    long a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
    long x, d, n = 1;
    int i, neg = 0;

    if(std::round(f*maxden)/maxden == 0) { den = 1; num = 0; return; }

    if(maxden <= 1) { den = 1; num = (long) f; return; }

    if(f < 0) { neg = 1; f = -f; }

    while (f != floor(f)) { n <<= 1; f *= 2; }
    d = f;

    /* continued fraction and check denominator each step */
    for (i = 0; i < 64; i++) {
        a = n ? d / n : 0;
        if (i && !a) break;

        x = d; d = n; n = x % n;

        x = a;
        if (k[1] * a + k[0] >= maxden) {
            x = (maxden - k[0]) / k[1];
            if (x * 2 >= a || k[1] >= maxden)
                i = 65;
            else
                break;
        }

        h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
        k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
    }
    den = k[1];
    num = neg ? -h[1] : h[1];
}

auto cleanRationalNumbers(double* vals, long size, long maxden) -> void
{
    for(long i = 0; i < size; ++i)
    {
        long num, den;
        fraction(vals[i], maxden, num, den);
        vals[i] = static_cast<double>(num)/den;
    }
}

} // namespace Reaktoro
