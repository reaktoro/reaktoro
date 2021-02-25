// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "MathUtils.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/QR>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

auto linearlyIndependentCols(MatrixXdConstRef A) -> Indices
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);
    Indices indices(I.data(), I.data() + rank);
    return indices;
}

auto linearlyIndependentRows(MatrixXdConstRef A) -> Indices
{
    const MatrixXd At = A.transpose();
    return linearlyIndependentCols(At);
}

auto linearlyIndependentCols(MatrixXdConstRef A, MatrixXdRef B) -> Indices
{
    Indices indices = linearlyIndependentCols(A);
    MatrixXd C(A.rows(), indices.size());
    for(unsigned i = 0; i < indices.size(); ++i)
        C.col(i) = A.col(indices[i]);
    B.noalias() = C;
    return indices;
}

auto linearlyIndependentRows(MatrixXdConstRef A, MatrixXdRef B) -> Indices
{
    Indices indices = linearlyIndependentRows(A);
    MatrixXd C(indices.size(), A.cols());
    for(unsigned i = 0; i < indices.size(); ++i)
        C.row(i) = A.row(indices[i]);
    B.noalias() = C;
    return indices;
}

auto inverseShermanMorrison(MatrixXdConstRef invA, VectorXdConstRef D) -> MatrixXd
{
    MatrixXd invM = invA;
    for(unsigned i = 0; i < D.rows(); ++i)
        invM = invM - (D[i]/(1 + D[i]*invM(i, i)))*invM.col(i)*invM.row(i);
    return invM;
}

/// Return the numerator and denominator of the rational number closest to `x`.
/// This methods expects `0 <= x <= 1`.
/// @param x The number for which the closest rational number is sought.
/// @param maxden The maximum denominator that the rational number can have.
auto farey(double x, unsigned maxden) -> std::tuple<long, long>
{
    long a = 0, b = 1;
    long c = 1, d = 1;
    while(b <= maxden && d <= maxden)
    {
        double mediant = double(a+c)/(b+d);
        if(x == mediant) {
            if(b + d <= maxden) return std::make_tuple(a+c, b+d);
            if(d > b) return std::make_tuple(c, d);
            return std::make_tuple(a, b);
        }
        if(x > mediant) {
            a = a+c;
            b = b+d;
        }
        else {
            c = a+c;
            d = b+d;
        }
    }

    return (b > maxden) ? std::make_tuple(c, d) : std::make_tuple(a, b);
}

auto rationalize(double x, unsigned maxden) -> std::tuple<long, long>
{
    long a, b, sign = (x >= 0) ? +1 : -1;
    if(std::abs(x) > 1.0) {
        std::tie(a, b) = farey(1.0/std::abs(x), maxden);
        return std::make_tuple(sign*b, a);
    }
    else {
        std::tie(a, b) = farey(std::abs(x), maxden);
        return std::make_tuple(sign*a, b);
    }
}

auto cleanRationalNumbers(double* vals, long size, long maxden) -> void
{
    long num, den;
    for(long i = 0; i < size; ++i)
    {
        std::tie(num, den) = rationalize(vals[i], maxden);
        vals[i] = static_cast<double>(num)/den;
    }
}

auto cleanRationalNumbers(MatrixXdRef A, long maxden) -> void
{
    cleanRationalNumbers(A.data(), A.size(), maxden);
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

auto dot3p(VectorXdConstRef x, VectorXdConstRef y, double s) -> double
{
    return dot3p_(x, y, s);
}

auto residual3p(MatrixXdConstRef A, VectorXdConstRef x, VectorXdConstRef b) -> VectorXd
{
    const auto m = A.rows();
    VectorXd r = zeros(m);
    for(int k = 0; k < m; ++k)
        r[k] = dot3p(A.row(k), x, -b[k]);
    return r;
}

} // namespace Reaktoro
