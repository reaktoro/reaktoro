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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

// Forward declaration
template<typename V, typename T, typename P, typename N>
class ChemicalVectorBase;

/// A type that represents a vector of chemical scalars and their derivatives.
/// A *chemical scalar* is a quantity that depends on temperature, pressure,
/// and molar amounts of species. A ChemicalScalar holds not only its value,
/// but also its temperature, pressure, and molar partial derivatives. A
/// ChemicalVector is a vector representation of a collection of ChemicalScalar
/// instances.
/// @see ThermoVector, ThermoScalar, ChemicalScalar
using ChemicalVector = ChemicalVectorBase<Vector,Vector,Vector,Matrix>;

/// A template base class to represent a vector of chemical scalars and their partial derivatives.
/// @see ThermoScalar, ThermoVector, ChemicalScalar, ChemicalVector
template<typename V, typename T, typename P, typename N>
class ChemicalVectorBase
{
public:
    /// The vector of chemical scalars
    V val;

    /// The vector of partial temperature derivatives of the chemical scalars
    T ddt;

    /// The vector of partial pressure derivatives of the chemical scalars
    P ddp;

    /// The matrix of partial molar derivatives of the chemical scalars
    N ddn;

    /// Return a ChemicalVector with zeros and zero derivatives.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    static auto Zero(Index nrows, Index nspecies) -> ChemicalVectorBase
    {
        return ChemicalVectorBase<V,T,P,N>(nrows, nspecies, 0.0);
    }

    /// Return a ChemicalVector with ones and zero derivatives.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    static auto One(Index nrows, Index nspecies) -> ChemicalVectorBase
    {
        return ChemicalVectorBase<V,T,P,N>(nrows, nspecies, 1.0);
    }

    /// Return a ChemicalVector with a given constant and zero derivatives.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    static auto Constant(Index nrows, Index nspecies, double val) -> ChemicalVectorBase
    {
        return ChemicalVectorBase<V,T,P,N>(nrows, nspecies, val);
    }

    /// Construct a default ChemicalVectorBase instance.
    ChemicalVectorBase() {}

    /// Construct a ChemicalVectorBase instance with given number of rows and species.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    ChemicalVectorBase(Index nrows, Index nspecies)
    : ChemicalVectorBase<V,T,P,N>(nrows, nspecies, 0.0) {}

    /// Construct a ChemicalVectorBase instance with number of rows equal the given number of species.
    /// @param nspecies The number of species for the molar derivatives
    explicit ChemicalVectorBase(Index nspecies)
    : ChemicalVectorBase<V,T,P,N>(nspecies, nspecies, 0.0) {}

    /// Construct a ChemicalVectorBase instance with given number of rows and species, and a constant value.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    /// @param val The constant value
    ChemicalVectorBase(Index nrows, Index nspecies, double val)
    : ChemicalVectorBase<V,T,P,N>(constants(nrows, val), zeros(nrows), zeros(nrows), zeros(nrows, nspecies)) {}

    /// Construct a ChemicalVectorBase instance with given values and derivatives.
    /// @param val The vector of values of the chemical scalars
    /// @param ddt The vector of partial temperature derivatives of the chemical scalars
    /// @param ddp The vector of partial pressure derivatives of the chemical scalars
    /// @param ddn The matrix of partial molar derivatives of the chemical scalars
    ChemicalVectorBase(const V& val, const T& ddt, const P& ddp, const N& ddn)
    : val(val), ddt(ddt), ddp(ddp), ddn(ddn) {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR, typename NR>
    ChemicalVectorBase(const ChemicalVectorBase<VR,TR,PR,NR>& other)
    : val(other.val), ddt(other.ddt), ddp(other.ddp), ddn(other.ddn) {}

    /// Return the number of rows in this ChemicalVectorBase instance.
    auto size() const -> Index
    {
        return val.size();
    }

    /// Resize this ChemicalVectorBase instance with new number of rows and number of species.
    /// @param nrows The new number of rows
    /// @param nspecies The new number of species for the molar derivatives
    auto resize(Index nrows, Index nspecies) -> void
    {
        val.resize(nrows);
        ddt.resize(nrows);
        ddp.resize(nrows);
        ddn.resize(nrows, nspecies);
    }

    /// Resize this ChemicalVectorBase instance with number of rows equal the number of species.
    /// @param nspecies The new number of species for the molar derivatives
    auto resize(Index nspecies) -> void
    {
        resize(nspecies, nspecies);
    }

    /// Assign another ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        val = other.val;
        ddt = other.ddt;
        ddp = other.ddp;
        ddn = other.ddn;
        return *this;
    }

    /// Assign a ChemicalScalarBase instance to this.
    template<typename VR, typename NR>
    auto operator=(const ChemicalScalarBase<VR,NR>& other) -> ChemicalVectorBase&
    {
        val.fill(other.val);
        ddt.fill(other.ddt);
        ddp.fill(other.ddp);
        for(auto i = 0; i < ddn.rows(); ++i) ddn.row(i) = other.ddn;
        return *this;
    }

    /// Assign a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.fill(other.val);
        ddt.fill(other.ddt);
        ddp.fill(other.ddp);
        ddn.fill(0.0);
        return *this;
    }

    /// Assign a scalar to this.
    auto operator=(double other) -> ChemicalVectorBase&
    {
        val.fill(other);
        ddt.fill(0.0);
        ddp.fill(0.0);
        ddn.fill(0.0);
        return *this;
    }

    /// Assign-addition of a ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator+=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        val += other.val;
        ddt += other.ddt;
        ddp += other.ddp;
        ddn += other.ddn;
        return *this;
    }

    /// Assign-addition of a ThermoVectorBase instance to this.
    template<typename VR, typename TR, typename PR>
    auto operator+=(const ThermoVectorBase<VR,TR,PR>& other) -> ChemicalVectorBase&
    {
        val += other.val;
        ddt += other.ddt;
        ddp += other.ddp;
        return *this;
    }

    /// Assign-addition of a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.array() += other.val;
        ddt.array() += other.ddt;
        ddp.array() += other.ddp;
        return *this;
    }

    /// Assign-addition of a scalar to this.
    auto operator+=(double other) -> ChemicalVectorBase&
    {
        val.array() += other;
        return *this;
    }

    /// Assign-subtraction of a ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator-=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        val -= other.val;
        ddt -= other.ddt;
        ddp -= other.ddp;
        ddn -= other.ddn;
        return *this;
    }

    /// Assign-subtraction of a ThermoVectorBase instance to this.
    template<typename VR, typename TR, typename PR>
    auto operator-=(const ThermoVectorBase<VR,TR,PR>& other) -> ChemicalVectorBase&
    {
        val -= other.val;
        ddt -= other.ddt;
        ddp -= other.ddp;
        return *this;
    }

    /// Assign-subtraction of a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.array() -= other.val;
        ddt.array() -= other.ddt;
        ddp.array() -= other.ddp;
        return *this;
    }

    /// Assign-subtraction of a scalar to this.
    auto operator-=(double other) -> ChemicalVectorBase&
    {
        val.array() -= other;
        return *this;
    }

    /// Assign-multiplication of a ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator*=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        ddt = diag(val) * other.ddt + diag(other.val) * ddt;
        ddp = diag(val) * other.ddp + diag(other.val) * ddp;
        ddn = diag(val) * other.ddn + diag(other.val) * ddn;
        val = diag(val) * other.val;
        return *this;
    }

    /// Assign-multiplication of a scalar to this.
    auto operator*=(double other) -> ChemicalVectorBase&
    {
        val *= other;
        ddt *= other;
        ddp *= other;
        ddn *= other;
        return *this;
    }

    /// Assign-division of a scalar to this.
    auto operator/=(double other) -> ChemicalVectorBase&
    {
        *this *= 1.0/other;
        return *this;
    }

    /// Return a ChemicalScalarBase with reference to the chemical scalar in a given row.
    auto operator[](Index irow) -> ChemicalScalarBase<double&, decltype(rowascol(ddn, irow))>
    {
        return {val[irow], ddt[irow], ddp[irow], rowascol(ddn, irow)};
    }

    /// Return a ChemicalScalarBase with const reference to the chemical scalar in a given row.
    auto operator[](Index irow) const -> ChemicalScalarBase<const double&, decltype(rowascol(ddn, irow))>
    {
        return {val[irow], ddt[irow], ddp[irow], rowascol(ddn, irow)};
    }

    /// Return a ChemicalScalarBase with reference to the chemical scalar in a given row.
    auto row(Index irow) -> ChemicalScalarBase<double&, decltype(rowascol(ddn, irow))>
    {
        return {val[irow], ddt[irow], ddp[irow], rowascol(ddn, irow)};
    }

    /// Return a ChemicalScalarBase with const reference to the chemical scalar in a given row.
    auto row(Index irow) const -> ChemicalScalarBase<const double&, decltype(rowascol(ddn, irow))>
    {
        return {val[irow], ddt[irow], ddp[irow], rowascol(ddn, irow)};
    }

    /// Return a reference of a row of this ChemicalVectorBase instance.
    auto row(Index irow, Index icol, Index ncols) -> ChemicalScalarBase<double&, decltype(ddn.row(irow).segment(icol, ncols))>
    {
        return {val[irow], ddt[irow], ddp[irow], ddn.row(irow).segment(icol, ncols)};
    }

    /// Return a const reference of a row of this ChemicalVectorBase instance.
    auto row(Index irow, Index icol, Index ncols) const -> ChemicalScalarBase<const double&, decltype(ddn.row(irow).segment(icol, ncols))>
    {
        return {val[irow], ddt[irow], ddp[irow], ddn.row(irow).segment(icol, ncols)};
    }

    /// Return a reference of a sequence of rows of this ChemicalVectorBase instance
    auto rows(Index irow, Index nrows) -> ChemicalVectorBase<decltype(Reaktoro::rows(val, irow, nrows)), decltype(Reaktoro::rows(ddt, irow, nrows)), decltype(Reaktoro::rows(ddp, irow, nrows)), decltype(Reaktoro::rows(ddn, irow, nrows))>
    {
        return {Reaktoro::rows(val, irow, nrows), Reaktoro::rows(ddt, irow, nrows), Reaktoro::rows(ddp, irow, nrows), Reaktoro::rows(ddn, irow, nrows)};
    }

    /// Return a const reference of a sequence of rows of this ChemicalVectorBase instance.
    auto rows(Index irow, Index nrows) const -> ChemicalVectorBase<decltype(Reaktoro::rows(val, irow, nrows)), decltype(Reaktoro::rows(ddt, irow, nrows)), decltype(Reaktoro::rows(ddp, irow, nrows)), decltype(Reaktoro::rows(ddn, irow, nrows))>
    {
        return {Reaktoro::rows(val, irow, nrows), Reaktoro::rows(ddt, irow, nrows), Reaktoro::rows(ddp, irow, nrows), Reaktoro::rows(ddn, irow, nrows)};
    }

    /// Return a reference of a sequence of rows of this ChemicalVectorBase instance.
    auto rows(Index irow, Index icol, Index nrows, Index ncols) -> ChemicalVectorBase<decltype(Reaktoro::rows(val, irow, nrows)), decltype(Reaktoro::rows(ddt, irow, nrows)), decltype(Reaktoro::rows(ddp, irow, nrows)), decltype(Reaktoro::block(ddn, irow, icol, nrows, ncols))>
    {
        return {Reaktoro::rows(val, irow, nrows), Reaktoro::rows(ddt, irow, nrows), Reaktoro::rows(ddp, irow, nrows), Reaktoro::block(ddn, irow, icol, nrows, ncols)};
    }

    /// Return a const reference of a sequence of rows of this ChemicalVectorBase instance.
    auto rows(Index irow, Index icol, Index nrows, Index ncols) const -> ChemicalVectorBase<decltype(Reaktoro::rows(val, irow, nrows)), decltype(Reaktoro::rows(ddt, irow, nrows)), decltype(Reaktoro::rows(ddp, irow, nrows)), decltype(Reaktoro::block(ddn, irow, icol, nrows, ncols))>
    {
        return {Reaktoro::rows(val, irow, nrows), Reaktoro::rows(ddt, irow, nrows), Reaktoro::rows(ddp, irow, nrows), Reaktoro::block(ddn, irow, icol, nrows, ncols)};
    }

    /// Return a reference of some rows of this ChemicalVectorBase instance.
    auto rows(const Indices& irows) -> ChemicalVector
    {
        return {Reaktoro::rows(val, irows), Reaktoro::rows(ddt, irows), Reaktoro::rows(ddp, irows), Reaktoro::rows(ddn, irows)};
    }

    /// Return a const reference of some rows of this ChemicalVectorBase instance.
    auto rows(const Indices& irows) const -> ChemicalVector
    {
        return {Reaktoro::rows(val, irows), Reaktoro::rows(ddt, irows), Reaktoro::rows(ddp, irows), Reaktoro::rows(ddn, irows)};
    }

    /// Return a reference of some rows and cols of this ChemicalVectorBase instance.
    auto rows(const Indices& irows, const Indices& icols) -> ChemicalVector
    {
        return {MatrixViewRows<V>(val, irows), MatrixViewRows<T>(ddt, irows), MatrixViewRows<P>(ddp, irows), MatrixViewRowsCols<N>(ddn, irows, icols)};
    }

    /// Return a const reference of some rows and cols of this ChemicalVectorBase instance.
    auto rows(const Indices& irows, const Indices& icols) const -> ChemicalVector
    {
        return {MatrixViewRowsConst<V>(val, irows), MatrixViewRowsConst<T>(ddt, irows), MatrixViewRowsConst<P>(ddp, irows), MatrixViewRowsColsConst<N>(ddn, irows, icols)};
    }
};

/// Return a ChemicalVector representation of a vector of molar composition of species.
/// @param n The molar composition vector of a collection of species.
template<typename Derived>
auto composition(const Eigen::MatrixBase<Derived>& n) -> ChemicalVectorBase<decltype(n), decltype(zeros(n.rows())), decltype(zeros(n.rows())), decltype(identity(n.rows(), n.rows()))>
{
    return {n, zeros(n.rows()), zeros(n.rows()), identity(n.rows(), n.rows())};
}

template<typename V, typename T, typename P, typename N>
auto operator+(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVectorBase<V,T,P,N>
{
    return l;
}

template<typename V, typename T, typename P, typename N>
auto operator-(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVectorBase<decltype(-l.val), decltype(-l.ddt), decltype(-l.ddp), decltype(-l.ddn)>
{
    return {-l.val, -l.ddt, -l.ddp, -l.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(l.val + r.val), decltype(l.ddt + r.ddt), decltype(l.ddp + r.ddp), decltype(l.ddn + r.ddn)>
{
    return {l.val + r.val, l.ddt + r.ddt, l.ddp + r.ddp, l.ddn + r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVectorBase<decltype(l.val + r.val), decltype(l.ddt + r.ddt), decltype(l.ddp + r.ddp), decltype(l.ddn)>
{
    return {l.val + r.val, l.ddt + r.ddt, l.ddp + r.ddp, l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val + r.val*ones(l.size())), decltype(l.ddt + r.ddt*ones(l.size())), decltype(l.ddp + r.ddp*ones(l.size())), decltype(l.ddn)>
{
    return {l.val + r.val*ones(l.size()), l.ddt + r.ddt*ones(l.size()), l.ddp + r.ddp*ones(l.size()), l.ddn};
}

template<typename VL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename V, typename T, typename P, typename N>
auto operator+(const ChemicalVectorBase<V,T,P,N>& l, const Vector& r) -> ChemicalVectorBase<decltype(l.val + r),T,P,N>
{
    return {l.val + r, l.ddt, l.ddp, l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator+(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(l.val - r.val), decltype(l.ddt - r.ddt), decltype(l.ddp - r.ddp), decltype(l.ddn - r.ddn)>
{
    return {l.val - r.val, l.ddt - r.ddt, l.ddp - r.ddp, l.ddn - r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVectorBase<decltype(l.val - r.val), decltype(l.ddt - r.ddt), decltype(l.ddp - r.ddp), decltype(l.ddn)>
{
    return {l.val - r.val, l.ddt - r.ddt, l.ddp - r.ddp, l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

template<typename VL, typename TL, typename PL, typename NL, typename VR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val - r.val*ones(l.size())), decltype(l.ddt - r.ddt*ones(l.size())), decltype(l.ddp - r.ddp*ones(l.size())), decltype(l.ddn)>
{
    return {l.val - r.val*ones(l.size()), l.ddt - r.ddt*ones(l.size()), l.ddp - r.ddp*ones(l.size()), l.ddn};
}

template<typename VL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return -(r - l);
}

template<typename V, typename T, typename P, typename N>
auto operator-(const ChemicalVectorBase<V,T,P,N>& l, const Vector& r) -> ChemicalVectorBase<decltype(l.val - r),T,P,N>
{
    return {l.val - r, l.ddt, l.ddp, l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator-(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

template<typename V, typename T, typename P, typename N>
auto operator*(double l, const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<decltype(l * r.val), decltype(l * r.ddt), decltype(l * r.ddp), decltype(l * r.ddn)>
{
    return {l * r.val, l * r.ddt, l * r.ddp, l * r.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator*(const ChemicalVectorBase<V,T,P,N>& l, double r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator*(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,T,P,N>& r) -> ChemicalVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddt + l.ddt * r.val), decltype(l.val * r.ddp + l.ddp * r.val), decltype(l.val * r.ddn)>
{
    return {l.val * r.val, l.val * r.ddt + l.ddt * r.val, l.val * r.ddp + l.ddp * r.val, l.val * r.ddn};
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator*(const ChemicalVectorBase<VL,T,P,N>& l, const ThermoScalarBase<VR>& r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename NR>
auto operator*(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalScalarBase<VR,NR>& r) -> ChemicalVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddt + l.ddt * r.val), decltype(l.val * r.ddp + l.ddp * r.val), decltype(l.val * tr(r.ddn) + l.ddn * r.val)>
{
    return {l.val * r.val, l.val * r.ddt + l.ddt * r.val, l.val * r.ddp + l.ddp * r.val, l.val * tr(r.ddn) + l.ddn * r.val};
}

template<typename VL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator*(const ChemicalScalarBase<VL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator%(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddt + diag(r.val) * l.ddt), decltype(diag(l.val) * r.ddp + diag(r.val) * l.ddp), decltype(diag(l.val) * r.ddn + diag(r.val) * l.ddn)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddt + diag(r.val) * l.ddt,
            diag(l.val) * r.ddp + diag(r.val) * l.ddp,
            diag(l.val) * r.ddn + diag(r.val) * l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddt + diag(r.val) * l.ddt), decltype(diag(l.val) * r.ddp + diag(r.val) * l.ddp), decltype(diag(l.val) * r.ddn)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddt + diag(r.val) * l.ddt,
            diag(l.val) * r.ddp + diag(r.val) * l.ddp,
            diag(l.val) * r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator%(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> decltype(r % l)
{
    return r % l;
}

template<typename V, typename T, typename P, typename N>
auto operator%(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<decltype(diag(l) * r.val), decltype(diag(l) * r.ddt), decltype(diag(l) * r.ddp), decltype(diag(l) * r.ddn)>
{
    return {diag(l) * r.val, diag(l) * r.ddt, diag(l) * r.ddp, diag(l) * r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL>
auto operator%(const ChemicalVectorBase<VL,TL,PL,NL>& l, const Vector& r) -> decltype(r % l)
{
    return r % l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator/(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVector
{
    const Vector tmp = 1.0/(r.val % r.val);
    return {l.val/r.val,
            diag(tmp) * (diag(r.val) * l.ddt - diag(l.val) * r.ddt),
            diag(tmp) * (diag(r.val) * l.ddp - diag(l.val) * r.ddp),
            diag(tmp) * (diag(r.val) * l.ddn - diag(l.val) * r.ddn)};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator/(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVector
{
    const Vector tmp = 1.0/(r.val % r.val);
    return {l.val/r.val,
            diag(tmp) * (diag(r.val) * l.ddt - diag(l.val) * r.ddt),
            diag(tmp) * (diag(r.val) * l.ddp - diag(l.val) * r.ddp),
            diag(tmp) * (diag(r.val) * l.ddn)};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename NR>
auto operator/(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalScalarBase<VR,NR>& r) -> ChemicalVector
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val,
            tmp * (r.val * l.ddt - l.val * r.ddt),
            tmp * (r.val * l.ddp - l.val * r.ddp),
            tmp * (r.val * l.ddn - l.val * tr(r.ddn))};
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator/(const ChemicalVectorBase<VL,T,P,N>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val/r.val), decltype(double() * (l.ddt * r.val - l.val * r.ddt)), decltype(double() * (l.ddp * r.val - l.val * r.ddp)), decltype(double() * (l.ddn * r.val))>
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val,
            tmp * (l.ddt * r.val - l.val * r.ddt),
            tmp * (l.ddp * r.val - l.val * r.ddp),
            tmp * (l.ddn * r.val)};
}

template<typename V, typename T, typename P, typename N>
auto operator/(const ChemicalVectorBase<V,T,P,N>& l, double r) -> decltype((1.0/r) * l)
{
    return (1.0/r) * l;
}

template<typename V, typename T, typename P, typename N>
auto sqrt(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp1 = sqrt(l.val);
    const Vector tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto pow(const ChemicalVectorBase<V,T,P,N>& l, double power) -> ChemicalVector
{
    const Vector tmp1 = pow(l.val, power);
    const Vector tmp2 = power * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto exp(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp = exp(l.val);
    return {tmp, diag(tmp) * l.ddt, diag(tmp) * l.ddp, diag(tmp) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto log(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp1 = log(l.val);
    const Vector tmp2 = 1.0/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto log10(const ChemicalVectorBase<V,T,P,N>& l) -> decltype(log(l)/double())
{
    const double ln10 = 2.302585092994046;
    return log(l)/ln10;
}

template<typename V, typename T, typename P, typename N>
auto sum(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalScalarBase<double, decltype(l.ddn.rowwise().sum())>
{
    return {l.val.sum(), l.ddt.sum(), l.ddp.sum(), l.ddn.rowwise().sum()};
}

template<typename V, typename T, typename P, typename N>
auto min(const ChemicalVectorBase<V,T,P,N>& l) -> double
{
    return min(l.val);
}

template<typename V, typename T, typename P, typename N>
auto max(const ChemicalVectorBase<V,T,P,N>& l) -> double
{
    return max(l.val);
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator<(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val < r.val;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator<=(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val <= r.val;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator>(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val > r.val;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator>=(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val >= r.val;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator==(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val == r.val;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator!=(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> bool
{
    return l.val == r.val;
}

inline auto operator<<(std::ostream& out, const ChemicalVector& vector) -> std::ostream&
{
    out << vector.val;
    return out;
}

} // namespace Reaktoro
