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

#include "MathUtils.hpp"

// Eigen includes
#include <Reaktoro/Eigen/QR>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

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

auto fraction(double fin, long maxden, long& num, long& den) -> void
{
    double f = fin;

    // Assert the given number is finite (neither NaN nor INF)
    Assert(std::isfinite(f), "Could not compute the rational fraction "
        "of given floating-point number.", "The given number `" +
        std::to_string(fin) + "` is not finite (i.e., it is either NaN or INF).");

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

        Assert(n != 0, "Could not compute the rational fraction "
            "of given floating-point number.", "The given number `" +
            std::to_string(fin) + "` must have been spoiled by round-off errors.");

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

auto cleanRationalNumbers(Matrix& A, long maxden) -> void
{
    cleanRationalNumbers(A.data(), A.size(), maxden);
}

auto cleanRationalNumbers(Vector& x, long maxden) -> void
{
    cleanRationalNumbers(x.data(), x.size(), maxden);
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
