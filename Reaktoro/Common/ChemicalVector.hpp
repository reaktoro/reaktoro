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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declaration
template<typename V, typename T, typename P, typename N>
class ChemicalVectorBase;

/// A type that represents a vector of chemical properties and their derivatives.
/// A *chemical scalar* is a quantity that depends on temperature, pressure,
/// and mole amounts of species. A ChemicalScalar holds not only its value,
/// but also its partial derivatives with respect to temperature, pressure,
/// and species amounts. A ChemicalVector is a vector representation of
/// a collection of ChemicalScalar instances.
/// @see ThermoVector, ThermoScalar, ChemicalScalar
using ChemicalVector = ChemicalVectorBase<Vector,Vector,Vector,Matrix>;

/// A type that represents a vector of chemical properties and their derivatives.
using ChemicalVectorRef = ChemicalVectorBase<VectorRef,VectorRef,VectorRef,MatrixRef>;

/// A type that represents a vector of chemical properties and their derivatives.
using ChemicalVectorConstRef = ChemicalVectorBase<VectorConstRef,VectorConstRef,VectorConstRef,MatrixConstRef>;

/// A template base class to represent a vector of chemical properties and their partial derivatives.
/// @see ThermoScalar, ThermoVector, ChemicalScalar, ChemicalVector
template<typename V, typename T, typename P, typename N>
class ChemicalVectorBase
{
public:
    /// The vector of chemical scalars
    V val;

    /// The vector of partial temperature derivatives of the chemical scalars
    T ddT;

    /// The vector of partial pressure derivatives of the chemical scalars
    P ddP;

    /// The matrix of partial mole derivatives of the chemical scalars
    N ddn;

    /// Construct a default ChemicalVectorBase instance.
    ChemicalVectorBase() {}

    /// Construct a ChemicalVectorBase instance with number of rows equal to given number of species.
    /// @param nspecies The number of nspecies in the chemical vector.
    explicit ChemicalVectorBase(Index nspecies)
    : ChemicalVectorBase(nspecies, nspecies) {}

    /// Construct a ChemicalVectorBase instance with given number of rows and species.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    ChemicalVectorBase(Index nrows, Index nspecies)
    : ChemicalVectorBase(zeros(nrows), zeros(nrows), zeros(nrows), zeros(nrows, nspecies)) {}

    /// Construct a ChemicalVectorBase instance with given number of rows and species, and a constant value.
    /// @param nrows The number of rows in the chemical vector
    /// @param nspecies The number of species for the molar derivatives
    /// @param val The constant value
    ChemicalVectorBase(Index nrows, Index nspecies, double val)
    : ChemicalVectorBase(constants(nrows, val), zeros(nrows), zeros(nrows), zeros(nrows, nspecies)) {}

    /// Construct a ChemicalVectorBase instance with given values and derivatives.
    /// @param val The vector of values of the chemical scalars
    /// @param ddT The vector of partial temperature derivatives of the chemical scalars
    /// @param ddP The vector of partial pressure derivatives of the chemical scalars
    /// @param ddn The matrix of partial mole derivatives of the chemical scalars
    ChemicalVectorBase(const V& val, const T& ddT, const P& ddP, const N& ddn)
    : val(val), ddT(ddT), ddP(ddP), ddn(ddn) {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR, typename NR>
    ChemicalVectorBase(ChemicalVectorBase<VR,TR,PR,NR>& other)
    : val(other.val), ddT(other.ddT), ddP(other.ddP), ddn(other.ddn) {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR, typename NR>
    ChemicalVectorBase(const ChemicalVectorBase<VR,TR,PR,NR>& other)
    : val(other.val), ddT(other.ddT), ddP(other.ddP), ddn(other.ddn) {}

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
        val = zeros(nrows);
        ddT = zeros(nrows);
        ddP = zeros(nrows);
        ddn = zeros(nrows, nspecies);
    }

    /// Resize this ChemicalVectorBase instance with number of rows equal the number of species.
    /// @param nspecies The new number of species for the molar derivatives
    auto resize(Index nspecies) -> void
    {
        resize(nspecies, nspecies);
    }

    /// Assign a ChemicalScalarBase instance to this.
    template<typename VR, typename NR>
    auto fill(const ChemicalScalarBase<VR, NR>& other) -> void
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
        for(auto i = 0; i < ddn.rows(); ++i) ddn.row(i) = other.ddn;
    }

    /// Assign a ThermoScalarBase instance to this.
    template<typename VR>
    auto fill(const ThermoScalarBase<VR>& other) -> void
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
        ddn.fill(0.0);
    }

    /// Assign a scalarsto this.
    auto fill(double value) -> void
	{
    	val.fill(value);
    	ddT.fill(0.0);
    	ddP.fill(0.0);
    	ddn.fill(0.0);
	}

    /// Assign another ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        val = other.val;
        ddT = other.ddT;
        ddP = other.ddP;
        ddn = other.ddn;
        return *this;
    }

    /// Assign a ChemicalScalarBase instance to this.
    template<typename VR, typename NR>
    auto operator=(const ChemicalScalarBase<VR,NR>& other) -> ChemicalVectorBase&
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
        for(auto i = 0; i < ddn.rows(); ++i) ddn.row(i) = other.ddn;
        return *this;
    }

    /// Assign a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
        ddn.fill(0.0);
        return *this;
    }

    /// Assign a scalar to this.
    auto operator=(double other) -> ChemicalVectorBase&
    {
        val.fill(other);
        ddT.fill(0.0);
        ddP.fill(0.0);
        ddn.fill(0.0);
        return *this;
    }

    /// Assign-addition of a ChemicalVectorBase instance to this.
    template<typename VR, typename TR, typename PR, typename NR>
    auto operator+=(const ChemicalVectorBase<VR,TR,PR,NR>& other) -> ChemicalVectorBase&
    {
        val += other.val;
        ddT += other.ddT;
        ddP += other.ddP;
        ddn += other.ddn;
        return *this;
    }

    /// Assign-addition of a ThermoVectorBase instance to this.
    template<typename VR, typename TR, typename PR>
    auto operator+=(const ThermoVectorBase<VR,TR,PR>& other) -> ChemicalVectorBase&
    {
        val += other.val;
        ddT += other.ddT;
        ddP += other.ddP;
        return *this;
    }

    /// Assign-addition of a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.array() += other.val;
        ddT.array() += other.ddT;
        ddP.array() += other.ddP;
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
        ddT -= other.ddT;
        ddP -= other.ddP;
        ddn -= other.ddn;
        return *this;
    }

    /// Assign-subtraction of a ThermoVectorBase instance to this.
    template<typename VR, typename TR, typename PR>
    auto operator-=(const ThermoVectorBase<VR,TR,PR>& other) -> ChemicalVectorBase&
    {
        val -= other.val;
        ddT -= other.ddT;
        ddP -= other.ddP;
        return *this;
    }

    /// Assign-subtraction of a ThermoScalarBase instance to this.
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ChemicalVectorBase&
    {
        val.array() -= other.val;
        ddT.array() -= other.ddT;
        ddP.array() -= other.ddP;
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
        ddT = diag(val) * other.ddT + diag(other.val) * ddT;
        ddP = diag(val) * other.ddP + diag(other.val) * ddP;
        ddn = diag(val) * other.ddn + diag(other.val) * ddn;
        val = diag(val) * other.val;
        return *this;
    }

    /// Assign-multiplication of a scalar to this.
    auto operator*=(double other) -> ChemicalVectorBase&
    {
        val *= other;
        ddT *= other;
        ddP *= other;
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
    auto operator[](Index irow) -> ChemicalScalarBase<double&, decltype(ddn.row(irow))>
    {
        return {val[irow], ddT[irow], ddP[irow], ddn.row(irow)};
    }

    /// Return a ChemicalScalarBase with const reference to the chemical scalar in a given row.
    auto operator[](Index irow) const -> ChemicalScalarBase<double, decltype(ddn.row(irow))>
    {
        return {val[irow], ddT[irow], ddP[irow], ddn.row(irow)};
    }

    /// Return a view of an interval of the ChemicalVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param nrows The number of rows in the view.
    auto view(Index irow, Index nrows) -> ChemicalVectorRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows), rowsmap(ddn, irow, nrows)};
    }

    /// Return a view of an interval of the ChemicalVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param nrows The number of rows in the view.
    auto view(Index irow, Index nrows) const -> ChemicalVectorConstRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows), rowsmap(ddn, irow, nrows)};
    }

    /// Return a view of an interval of the ChemicalVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param icol The index of the col starting the view of the ddn derivatives.
    /// @param nrows The number of rows in the view.
    /// @param ncols The number of columns of the view of the ddn derivatives.
    auto view(Index irow, Index icol, Index nrows, Index ncols) -> ChemicalVectorRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows), blockmap(ddn, irow, icol, nrows, ncols)};
    }

    /// Return a view of an interval of the ChemicalVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param icol The index of the col starting the view of the ddn derivatives.
    /// @param nrows The number of rows in the view.
    /// @param ncols The number of columns of the view of the ddn derivatives.
    auto view(Index irow, Index icol, Index nrows, Index ncols) const -> ChemicalVectorConstRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows), blockmap(ddn, irow, icol, nrows, ncols)};
    }

    /// Explicitly converts this ChemicalVector instance into a Vector.
    explicit operator Vector() const
    {
        return val;
    }
};

/// Return a ChemicalVectorBase expression representing zeros with same dimension of given vector.
template<typename V, typename T, typename P, typename N>
auto zeros(const ChemicalVectorBase<V,T,P,N>& v) -> ChemicalVectorBase<decltype(zeros(0)), decltype(zeros(0)), decltype(zeros(0)), decltype(zeros(0,0))>
{
    const Index n = v.size();
    return {zeros(n), zeros(n), zeros(n), zeros(n, n)};
}

/// Return a ChemicalVectorBase expression representing ones with same dimension of given vector.
template<typename V, typename T, typename P, typename N>
auto ones(const ChemicalVectorBase<V,T,P,N>& v) -> ChemicalVectorBase<decltype(ones(0)), decltype(zeros(0)), decltype(zeros(0)), decltype(zeros(0,0))>
{
    const Index n = v.size();
    return {ones(n), zeros(n), zeros(n), zeros(n, n)};
}

/// A type that describes temperature in units of K
class Composition : public ChemicalVectorBase<VectorConstRef, decltype(zeros(0)), decltype(zeros(0)), decltype(identity(0,0))>
{
public:
	/// Auxiliary base type
	using Base = ChemicalVectorBase<VectorConstRef, decltype(zeros(0)), decltype(zeros(0)), decltype(identity(0,0))>;

    /// Construct a Composition instance with given composition vector.
    Composition(VectorConstRef n) : Base(n, zeros(n.rows()), zeros(n.rows()), identity(n.rows(), n.rows())) {}
};

template<typename V, typename T, typename P, typename N>
auto operator<<(std::ostream& out, const ChemicalVectorBase<V,T,P,N>& r) -> std::ostream&
{
    out << r.val;
    return out;
}

template<typename V, typename T, typename P, typename N>
auto operator+(const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<V,T,P,N>
{
    return r;
}

template<typename V, typename T, typename P, typename N>
auto operator-(const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<decltype(-r.val), decltype(-r.ddT), decltype(-r.ddP), decltype(-r.ddn)>
{
    return {-r.val, -r.ddT, -r.ddP, -r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(l.val + r.val), decltype(l.ddT + r.ddT), decltype(l.ddP + r.ddP), decltype(l.ddn + r.ddn)>
{
    return {l.val + r.val, l.ddT + r.ddT, l.ddP + r.ddP, l.ddn + r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVectorBase<decltype(l.val + r.val), decltype(l.ddT + r.ddT), decltype(l.ddP + r.ddP), decltype(l.ddn)>
{
    return {l.val + r.val, l.ddT + r.ddT, l.ddP + r.ddP, l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR>
auto operator+(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val + r.val*ones(l.size())), decltype(l.ddT + r.ddT*ones(l.size())), decltype(l.ddP + r.ddP*ones(l.size())), decltype(l.ddn)>
{
    return {l.val + r.val*ones(l.size()), l.ddT + r.ddT*ones(l.size()), l.ddP + r.ddP*ones(l.size()), l.ddn};
}

template<typename VL, typename VR, typename TR, typename PR, typename NR>
auto operator+(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename V, typename T, typename P, typename N>
auto operator+(const ChemicalVectorBase<V,T,P,N>& l, const Vector& r) -> ChemicalVectorBase<decltype(l.val + r),T,P,N>
{
    return {l.val + r, l.ddT, l.ddP, l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator+(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(l.val - r.val), decltype(l.ddT - r.ddT), decltype(l.ddP - r.ddP), decltype(l.ddn - r.ddn)>
{
    return {l.val - r.val, l.ddT - r.ddT, l.ddP - r.ddP, l.ddn - r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVectorBase<decltype(l.val - r.val), decltype(l.ddT - r.ddT), decltype(l.ddP - r.ddP), decltype(l.ddn)>
{
    return {l.val - r.val, l.ddT - r.ddT, l.ddP - r.ddP, l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

template<typename VL, typename TL, typename PL, typename NL, typename VR>
auto operator-(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val - r.val*ones(l.size())), decltype(l.ddT - r.ddT*ones(l.size())), decltype(l.ddP - r.ddP*ones(l.size())), decltype(l.ddn)>
{
    return {l.val - r.val*ones(l.size()), l.ddT - r.ddT*ones(l.size()), l.ddP - r.ddP*ones(l.size()), l.ddn};
}

template<typename VL, typename VR, typename TR, typename PR, typename NR>
auto operator-(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r + l)
{
    return -(r - l);
}

template<typename V, typename T, typename P, typename N>
auto operator-(const ChemicalVectorBase<V,T,P,N>& l, const Vector& r) -> ChemicalVectorBase<decltype(l.val - r),T,P,N>
{
    return {l.val - r, l.ddT, l.ddP, l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator-(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

template<typename V, typename T, typename P, typename N>
auto operator*(double l, const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<decltype(l * r.val), decltype(l * r.ddT), decltype(l * r.ddP), decltype(l * r.ddn)>
{
    return {l * r.val, l * r.ddT, l * r.ddP, l * r.ddn};
}

template<typename V, typename T, typename P, typename N>
auto operator*(const ChemicalVectorBase<V,T,P,N>& l, double r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator*(const ThermoScalarBase<VL>& l, const ChemicalVectorBase<VR,T,P,N>& r) -> ChemicalVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddT + l.ddT * r.val), decltype(l.val * r.ddP + l.ddP * r.val), decltype(l.val * r.ddn)>
{
    return {l.val * r.val, l.val * r.ddT + l.ddT * r.val, l.val * r.ddP + l.ddP * r.val, l.val * r.ddn};
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator*(const ChemicalVectorBase<VL,T,P,N>& l, const ThermoScalarBase<VR>& r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename NR>
auto operator*(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalScalarBase<VR,NR>& r) -> ChemicalVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddT + l.ddT * r.val), decltype(l.val * r.ddP + l.ddP * r.val), decltype(l.val * r.ddn + l.ddn * r.val)>
{
    return {l.val * r.val, l.val * r.ddT + l.ddT * r.val, l.val * r.ddP + l.ddP * r.val, l.val * r.ddn + l.ddn * r.val};
}

template<typename VL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator*(const ChemicalScalarBase<VL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator%(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddT + diag(r.val) * l.ddT), decltype(diag(l.val) * r.ddP + diag(r.val) * l.ddP), decltype(diag(l.val) * r.ddn + diag(r.val) * l.ddn)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddT + diag(r.val) * l.ddT,
            diag(l.val) * r.ddP + diag(r.val) * l.ddP,
            diag(l.val) * r.ddn + diag(r.val) * l.ddn};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR, typename NR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddT + diag(r.val) * l.ddT), decltype(diag(l.val) * r.ddP + diag(r.val) * l.ddP), decltype(diag(l.val) * r.ddn)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddT + diag(r.val) * l.ddT,
            diag(l.val) * r.ddP + diag(r.val) * l.ddP,
            diag(l.val) * r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR>
auto operator%(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> decltype(r % l)
{
    return r % l;
}

template<typename V, typename T, typename P, typename N>
auto operator%(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalVectorBase<decltype(diag(l) * r.val), decltype(diag(l) * r.ddT), decltype(diag(l) * r.ddP), decltype(diag(l) * r.ddn)>
{
    return {diag(l) * r.val, diag(l) * r.ddT, diag(l) * r.ddP, diag(l) * r.ddn};
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
            diag(tmp) * (diag(r.val) * l.ddT - diag(l.val) * r.ddT),
            diag(tmp) * (diag(r.val) * l.ddP - diag(l.val) * r.ddP),
            diag(tmp) * (diag(r.val) * l.ddn - diag(l.val) * r.ddn)};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto operator/(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ChemicalVector
{
    const Vector tmp = 1.0/(r.val % r.val);
    return {l.val/r.val,
            diag(tmp) * (diag(r.val) * l.ddT - diag(l.val) * r.ddT),
            diag(tmp) * (diag(r.val) * l.ddP - diag(l.val) * r.ddP),
            diag(tmp) * (diag(r.val) * l.ddn)};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename NR>
auto operator/(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalScalarBase<VR,NR>& r) -> ChemicalVector
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val,
            tmp * (r.val * l.ddT - l.val * r.ddT),
            tmp * (r.val * l.ddP - l.val * r.ddP),
            tmp * (r.val * l.ddn - l.val * r.ddn)};
}

template<typename VL, typename VR, typename T, typename P, typename N>
auto operator/(const ChemicalVectorBase<VL,T,P,N>& l, const ThermoScalarBase<VR>& r) -> ChemicalVectorBase<decltype(l.val/r.val), decltype(double() * (l.ddT * r.val - l.val * r.ddT)), decltype(double() * (l.ddP * r.val - l.val * r.ddP)), decltype(double() * (l.ddn * r.val))>
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val,
            tmp * (l.ddT * r.val - l.val * r.ddT),
            tmp * (l.ddP * r.val - l.val * r.ddP),
            tmp * (l.ddn * r.val)};
}

template<typename V, typename T, typename P, typename N>
auto operator/(const ChemicalVectorBase<V,T,P,N>& l, double r) -> decltype((1.0/r) * l)
{
    return (1.0/r) * l;
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

/// Return a ChemicalScalarBase with reference to the chemical scalar in a given row.
template<typename V, typename T, typename P, typename N>
auto row(ChemicalVectorBase<V,T,P,N>& vec, Index irow) -> ChemicalScalarBase<double&, decltype(vec.ddn.row(irow))>
{
	return {vec.val[irow], vec.ddT[irow], vec.ddP[irow], vec.ddn.row(irow)};
}

/// Return a ChemicalScalarBase with const reference to the chemical scalar in a given row.
template<typename V, typename T, typename P, typename N>
auto row(const ChemicalVectorBase<V,T,P,N>& vec, Index irow) -> ChemicalScalarBase<const double&, decltype(vec.ddn.row(irow))>
{
	return {vec.val[irow], vec.ddT[irow], vec.ddP[irow], vec.ddn.row(irow)};
}

/// Return a reference of a row of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto row(ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index icol, Index ncols) -> ChemicalScalarBase<double&, decltype(vec.ddn.row(irow).segment(icol, ncols))>
{
	return {vec.val[irow], vec.ddT[irow], vec.ddP[irow], vec.ddn.row(irow).segment(icol, ncols)};
}

/// Return a const reference of a row of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto row(const ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index icol, Index ncols) -> ChemicalScalarBase<const double&, decltype(vec.ddn.row(irow).segment(icol, ncols))>
{
	return {vec.val[irow], vec.ddT[irow], vec.ddP[irow], vec.ddn.row(irow).segment(icol, ncols)};
}

/// Return a reference of a sequence of rows of this ChemicalVectorBase instance
template<typename V, typename T, typename P, typename N>
auto rows(ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index nrows) -> ChemicalVectorBase<decltype(rows(vec.val, irow, nrows)), decltype(rows(vec.ddT, irow, nrows)), decltype(rows(vec.ddP, irow, nrows)), decltype(rows(vec.ddn, irow, nrows))>
{
	return {rows(vec.val, irow, nrows), rows(vec.ddT, irow, nrows), rows(vec.ddP, irow, nrows), rows(vec.ddn, irow, nrows)};
}

/// Return a const reference of a sequence of rows of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(const ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index nrows) -> ChemicalVectorBase<decltype(rows(vec.val, irow, nrows)), decltype(rows(vec.ddT, irow, nrows)), decltype(rows(vec.ddP, irow, nrows)), decltype(rows(vec.ddn, irow, nrows))>
{
	return {rows(vec.val, irow, nrows), rows(vec.ddT, irow, nrows), rows(vec.ddP, irow, nrows), rows(vec.ddn, irow, nrows)};
}

/// Return a reference of a sequence of rows of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index icol, Index nrows, Index ncols) -> ChemicalVectorBase<decltype(rows(vec.val, irow, nrows)), decltype(rows(vec.ddT, irow, nrows)), decltype(rows(vec.ddP, irow, nrows)), decltype(block(vec.ddn, irow, icol, nrows, ncols))>
{
	return {rows(vec.val, irow, nrows), rows(vec.ddT, irow, nrows), rows(vec.ddP, irow, nrows), block(vec.ddn, irow, icol, nrows, ncols)};
}

/// Return a const reference of a sequence of rows of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(const ChemicalVectorBase<V,T,P,N>& vec, Index irow, Index icol, Index nrows, Index ncols) -> ChemicalVectorBase<decltype(rows(vec.val, irow, nrows)), decltype(rows(vec.ddT, irow, nrows)), decltype(rows(vec.ddP, irow, nrows)), decltype(block(vec.ddn, irow, icol, nrows, ncols))>
{
	return {rows(vec.val, irow, nrows), rows(vec.ddT, irow, nrows), rows(vec.ddP, irow, nrows), block(vec.ddn, irow, icol, nrows, ncols)};
}

/// Return a reference of some rows of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(ChemicalVectorBase<V,T,P,N>& vec, const Indices& irows) -> ChemicalVector
{
	return {rows(vec.val, irows), rows(vec.ddT, irows), rows(vec.ddP, irows), rows(vec.ddn, irows)};
}

/// Return a const reference of some rows of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(const ChemicalVectorBase<V,T,P,N>& vec, const Indices& irows) -> ChemicalVector
{
	return {rows(vec.val, irows), rows(vec.ddT, irows), rows(vec.ddP, irows), rows(vec.ddn, irows)};
}

/// Return a reference of some rows and cols of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(ChemicalVectorBase<V,T,P,N>& vec, const Indices& irows, const Indices& icols) -> ChemicalVector
{
	return {rows(vec.val, irows), rows(vec.ddT, irows), rows(vec.ddP, irows), submatrix(vec.ddn, irows, icols)};
}

/// Return a const reference of some rows and cols of this ChemicalVectorBase instance.
template<typename V, typename T, typename P, typename N>
auto rows(const ChemicalVectorBase<V,T,P,N>& vec, const Indices& irows, const Indices& icols) -> ChemicalVector
{
	return {rows(vec.val, irows), rows(vec.ddT, irows), rows(vec.ddP, irows), submatrix(vec.ddn, irows, icols)};
}

template<typename V, typename T, typename P, typename N>
auto abs(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp1 = abs(l.val);
    const Vector tmp2 = l.val/tmp1;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto sqrt(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp1 = sqrt(l.val);
    const Vector tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto pow(const ChemicalVectorBase<V,T,P,N>& l, double power) -> ChemicalVector
{
    const Vector tmp1 = pow(l.val, power);
    const Vector tmp2 = power * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto exp(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp = exp(l.val);
    return {tmp, diag(tmp) * l.ddT, diag(tmp) * l.ddP, diag(tmp) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto log(const ChemicalVectorBase<V,T,P,N>& l) -> ChemicalVector
{
    const Vector tmp1 = log(l.val);
    const Vector tmp2 = 1.0/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP, diag(tmp2) * l.ddn};
}

template<typename V, typename T, typename P, typename N>
auto log10(const ChemicalVectorBase<V,T,P,N>& l) -> decltype(log(l)/double())
{
    const double ln10 = 2.302585092994046;
    return log(l)/ln10;
}

template<typename V, typename T, typename P, typename N>
auto min(const ChemicalVectorBase<V,T,P,N>& r) -> double
{
    return min(r.val);
}

template<typename V, typename T, typename P, typename N>
auto max(const ChemicalVectorBase<V,T,P,N>& r) -> double
{
    return max(r.val);
}

template<typename V, typename T, typename P, typename N>
auto sum(const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalScalarBase<double, decltype(r.ddn.colwise().sum())>
{
    return {r.val.sum(), r.ddT.sum(), r.ddP.sum(), r.ddn.colwise().sum()};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto dot(const ChemicalVectorBase<VL,TL,PL,NL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalScalarBase<double, decltype(tr(l.val) * r.ddn + tr(r.val) * l.ddn)>
{
    return {tr(l.val) * r.val, tr(l.val) * r.ddT + tr(r.val) * l.ddT, tr(l.val) * r.ddP + tr(r.val) * l.ddP, tr(l.val) * r.ddn + tr(r.val) * l.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto dot(const ThermoVectorBase<VL,TL,PL>& l, const ChemicalVectorBase<VR,TR,PR,NR>& r) -> ChemicalScalarBase<double, decltype(tr(l.val) * r.ddn)>
{
    return {tr(l.val) * r.val, tr(l.val) * r.ddT + tr(r.val) * l.ddT, tr(l.val) * r.ddP + tr(r.val) * l.ddP, tr(l.val) * r.ddn};
}

template<typename VL, typename TL, typename PL, typename NL, typename VR, typename TR, typename PR, typename NR>
auto dot(const ChemicalVectorBase<VR,TR,PR,NR>& l, const ThermoVectorBase<VL,TL,PL>& r) -> decltype(dot(r, l))
{
    return dot(r, l);
}

template<typename V, typename T, typename P, typename N>
auto dot(const Vector& l, const ChemicalVectorBase<V,T,P,N>& r) -> ChemicalScalarBase<double, decltype(tr(l) * r.ddn)>
{
    return {tr(l) * r.val, tr(l) * r.ddT, tr(l) * r.ddP, tr(l) * r.ddn};
}

template<typename V, typename T, typename P, typename N>
auto dot(const ChemicalVectorBase<V,T,P,N>& l, const Vector& r) -> decltype(dot(r, l))
{
    return dot(r, l);
}

} // namespace Reaktoro
