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
    lu(A, L, U, P, Q);
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
auto luSquare(const MatrixType& A, Matrix& L, Matrix& U, PermutationMatrix& P, PermutationMatrix& Q) -> void
{
    const Index m = A.rows();
    Eigen::PartialPivLU<Matrix> lu = A.lu();
    L = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitLower>();
    U = lu.matrixLU().triangularView<Eigen::Upper>();
    P = lu.permutationP();
    Q.setIdentity(m);
}

/// Compute LU decomposition of a rectangular matrix (more rows than cols).
template<typename MatrixType>
auto luRectangular1(const MatrixType& A, Matrix& L, Matrix& U, PermutationMatrix& P, PermutationMatrix& Q) -> void
{
    const Index m = A.rows();
    const Index n = A.cols();
    Matrix M(n, n);
    M.topRows(m) = A;
    M.bottomLeftCorner(n-m, m).fill(0.0);
    M.bottomRightCorner(n-m, n-m).setIdentity();
    Eigen::PartialPivLU<Matrix> lu = M.lu();
    L = lu.matrixLU().topLeftCorner(m, m).triangularView<Eigen::UnitLower>();
    U = lu.matrixLU().topRows(m).triangularView<Eigen::Upper>();
    P.resize(m);
    std::copy(lu.permutationP().indices().data(),
              lu.permutationP().indices().data() + m,
              P.indices().data());
    Q.setIdentity(n);
    for(Index i = 0; i < m; ++i)
        for(Index j = i; j < n; ++j)
            if(U(i, j) != 0.0)
            {
                auto tmp = Q.indices()[i];
                Q.indices()[i] = Q.indices()[j];
                Q.indices()[j] = tmp;
                break;
            }
    U = U * Q;
//    std::cout << "Q = \n" << Q.indices() << std::endl;
}

/// Compute LU decomposition of a rectangular matrix (more cols than rows).
template<typename MatrixType>
auto luRectangular2(const MatrixType& A, Matrix& L, Matrix& U, PermutationMatrix& P, PermutationMatrix& Q) -> void
{
    const Index m = A.rows();
    const Index n = A.cols();
    Matrix M(m, m);
    M.leftCols(n) = A;
    M.topRightCorner(n, m-n).fill(0.0);
    M.bottomRightCorner(m-n, m-n).setIdentity();
    Eigen::PartialPivLU<Matrix> lu = M.lu();
    L = lu.matrixLU().leftCols(n).triangularView<Eigen::UnitLower>();
    U = lu.matrixLU().topLeftCorner(n, n).triangularView<Eigen::Upper>();
    P = lu.permutationP();
}

auto lu(const Matrix& A, Matrix& L, Matrix& U, PermutationMatrix& P, PermutationMatrix& Q) -> void
{
    const Index m = A.rows();
    const Index n = A.cols();
    if(m == n) luSquare(A, L, U, P, Q);
    else if(m < n) luRectangular1(A, L, U, P, Q);
    else luRectangular2(A, L, U, P, Q);
}

template<typename VectorTypeX, typename VectorTypeY>
auto dot3p_(const VectorTypeX& x, const VectorTypeY& y, double s) -> double
{
   double shi = double(float(s));
   double slo = s - shi;
   for(int k = 0; k < x.size(); ++k)
   {
      double xhi = double(float(x[k]));
      double xlo = x[k] - xhi;
      double yhi = double(float(y[k]));
      double ylo = y[k] - yhi;
      double tmp = xhi*yhi;
      double zhi = double(float(tmp));
      double zlo = tmp - zhi + xhi*ylo + xlo*yhi + xlo*ylo;

      tmp = shi + zhi;
      double del = tmp - shi - zhi;
      shi = double(float(tmp));
      slo = tmp - shi + slo + zlo - del;
   }

   s = shi + slo;
   return s;
}

auto dot3p(const Vector& x, const Vector& y, double s) -> double
{
    return dot3p_(x, y, s);
}

auto residual3p(const Matrix& A, const Vector& x, const Vector& b) -> Vector
{
    const auto m = A.rows();
    Vector r = zeros(m);
    for(int k = 0; k < m; ++k)
        r[k] = dot3p(A.row(k), x, -b[k]);
    return r;
}

} // namespace Reaktoro
