// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#pragma once

#define EIGEN_MATRIX_PLUGIN "Reaktoro/Math/EigenMatrixPlugin.hpp"

// Eigen includes
#include <Reaktoro/Math/Eigen/Core>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Eigen {

template<typename Derived, typename Indices>
class MatrixRowsView;

template<typename Derived, typename Indices>
class MatrixRowsViewConst;

template<typename Derived, typename Indices>
class MatrixColsView;

template<typename Derived, typename Indices>
class MatrixColsViewConst;

template<typename Derived, typename Indices>
class MatrixSubView;

template<typename Derived, typename Indices>
class MatrixSubViewConst;

namespace internal {

template<typename Derived, typename Indices>
struct traits<MatrixRowsView<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

template<typename Derived, typename Indices>
struct traits<MatrixRowsViewConst<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

template<typename Derived, typename Indices>
struct traits<MatrixColsView<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

template<typename Derived, typename Indices>
struct traits<MatrixColsViewConst<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

template<typename Derived, typename Indices>
struct traits<MatrixSubView<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

template<typename Derived, typename Indices>
struct traits<MatrixSubViewConst<Derived, Indices>>
{
	typedef Eigen::Dense StorageKind;
	typedef Eigen::MatrixXpr XprKind;
	typedef typename Derived::Scalar Scalar;
	typedef typename Derived::Index Index;
	enum {
		Flags = Eigen::ColMajor | EvalBeforeNestingBit | EvalBeforeAssigningBit,
		RowsAtCompileTime = Derived::RowsAtCompileTime,
		ColsAtCompileTime = Derived::ColsAtCompileTime,
		MaxRowsAtCompileTime = Derived::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = Derived::MaxColsAtCompileTime,
		CoeffReadCost = Derived::CoeffReadCost
	};
};

} // namespace internal

template<typename Derived, typename Indices>
class MatrixRowsView : public MatrixBase<MatrixRowsView<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixRowsView> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixRowsView)

    MatrixRowsView(MatrixBase<Derived>& mat, const Indices& irows)
      : m_mat(mat), m_irows(irows)
    {}

    template<typename DerivedOther>
    auto operator=(const MatrixBase<DerivedOther>& other) -> MatrixRowsView&
	{
    	for(Index i = 0; i < rows(); ++i)
			for(Index j = 0; j < cols(); ++j)
				coeff(i, j) = other(i, j);
    	return *this;
	}

    auto rows() const -> Index { return m_irows.size(); }
    auto cols() const -> Index { return m_mat.cols(); }

    auto coeff(Index row, Index col) -> Scalar& { return m_mat(m_irows[row], col); }
    auto coeff(Index row, Index col) const -> Scalar { return m_mat(m_irows[row], col); }
    auto operator()(Index row, Index col) -> Scalar& { return coeff(row, col); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    MatrixBase<Derived>& m_mat;
    Indices m_irows;
};

template<typename Derived, typename Indices>
class MatrixRowsViewConst : public MatrixBase<MatrixRowsViewConst<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixRowsViewConst> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixRowsViewConst)

    MatrixRowsViewConst(const MatrixBase<Derived>& mat, const Indices& irows)
      : m_mat(mat), m_irows(irows)
    {}

    auto rows() const -> Index { return m_irows.size(); }
    auto cols() const -> Index { return m_mat.cols(); }

    auto coeff(Index row, Index col) const -> Scalar { return m_mat(m_irows[row], col); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    const MatrixBase<Derived>& m_mat;
    Indices m_irows;
};

template<typename Derived, typename Indices>
class MatrixColsView : public MatrixBase<MatrixColsView<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixColsView> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixColsView)

    MatrixColsView(MatrixBase<Derived>& mat, const Indices& icols)
      : m_mat(mat), m_icols(icols)
    {}

    template<typename DerivedOther>
    auto operator=(const MatrixBase<DerivedOther>& other) -> MatrixColsView&
	{
    	for(Index i = 0; i < rows(); ++i)
			for(Index j = 0; j < cols(); ++j)
				coeff(i, j) = other(i, j);
    	return *this;
	}

    auto rows() const -> Index { return m_mat.rows(); }
    auto cols() const -> Index { return m_icols.size(); }

    auto coeff(Index row, Index col) -> Scalar& { return m_mat(row, m_icols[col]); }
    auto coeff(Index row, Index col) const -> Scalar { return m_mat(row, m_icols[col]); }
    auto operator()(Index row, Index col) -> Scalar& { return coeff(row, col); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    MatrixBase<Derived>& m_mat;
    Indices m_icols;
};

template<typename Derived, typename Indices>
class MatrixColsViewConst : public MatrixBase<MatrixColsViewConst<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixColsViewConst> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixColsViewConst)

    MatrixColsViewConst(const MatrixBase<Derived>& mat, const Indices& icols)
      : m_mat(mat), m_icols(icols)
    {}

    auto rows() const -> Index { return m_mat.rows(); }
    auto cols() const -> Index { return m_icols.size(); }

    auto coeff(Index row, Index col) const -> Scalar { return m_mat(row, m_icols[col]); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    const MatrixBase<Derived>& m_mat;
    Indices m_icols;
};

template<typename Derived, typename Indices>
class MatrixSubView : public MatrixBase<MatrixSubView<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixSubView> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixSubView)

    MatrixSubView(MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols)
      : m_mat(mat), m_irows(irows), m_icols(icols)
    {}

    template<typename DerivedOther>
    auto operator=(const MatrixBase<DerivedOther>& other) -> MatrixSubView&
	{
    	for(Index i = 0; i < rows(); ++i)
			for(Index j = 0; j < cols(); ++j)
				coeff(i, j) = other(i, j);
    	return *this;
	}

    auto rows() const -> Index { return m_irows.size(); }
    auto cols() const -> Index { return m_icols.size(); }

    auto coeff(Index row, Index col) -> Scalar& { return m_mat(m_irows[row], m_icols[col]); }
    auto coeff(Index row, Index col) const -> Scalar { return m_mat(m_irows[row], m_icols[col]); }
    auto operator()(Index row, Index col) -> Scalar& { return coeff(row, col); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    MatrixBase<Derived>& m_mat;
    Indices m_irows, m_icols;
};

template<typename Derived, typename Indices>
class MatrixSubViewConst : public MatrixBase<MatrixSubViewConst<Derived, Indices>>
{
public:
    typedef MatrixBase<MatrixSubViewConst> Base;
    typedef typename Derived::PlainObject PlainObject;
    EIGEN_DENSE_PUBLIC_INTERFACE(MatrixSubViewConst)

    MatrixSubViewConst(const MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols)
      : m_mat(mat), m_irows(irows), m_icols(icols)
    {}

    auto rows() const -> Index { return m_irows.size(); }
    auto cols() const -> Index { return m_icols.size(); }

    auto coeff(Index row, Index col) const -> Scalar { return m_mat(m_irows[row], m_icols[col]); }
    auto operator()(Index row, Index col) const -> Scalar { return coeff(row, col); }

    operator PlainObject() const
    {
    	PlainObject res(rows(), cols());
    	for(Index i = 0; i < rows(); ++i)
    		for(Index j = 0; j < cols(); ++j)
    			res(i, j) = coeff(i, j);
    	return res;
    }

private:
    const MatrixBase<Derived>& m_mat;
    Indices m_irows, m_icols;
};

} // namespace Eigen

namespace Reaktoro {

/// Define an alias to the vector type of the Eigen library
using Vector = Eigen::VectorXd;

/// Define an alias to the matrix type of the Eigen library
using Matrix = Eigen::MatrixXd;

/// Define an alias to a permutation matrix type of the Eigen library
using PermutationMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
auto zeros(Index rows) -> decltype(Vector::Zero(rows));

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
auto ones(Index rows) -> decltype(Vector::Ones(rows));

/// Return an expression of a vector with random entries
/// @param rows The number of rows
/// @return The expression of a vector with random entries equal to one
auto random(Index rows) -> decltype(Vector::Random(rows));

/// Return a linearly spaced vector
/// @param rows The number of rows
/// @param start The start of the sequence
/// @param stop The stop of the sequence
/// @return The expression of a vector with linearly spaced entries
auto linspace(Index rows, double start, double stop) -> decltype(Vector::LinSpaced(rows, start, stop));

/// Return an expression of a unit vector
/// @param rows The number of rows
/// @param i The index at which the component is one
/// @return The expression of a unit vector
auto unit(Index rows, Index i) -> decltype(Vector::Unit(rows, i));

/// Return an expression of a zero matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a zero matrix
auto zeros(Index rows, Index cols) -> decltype(Matrix::Zero(rows, cols));

/// Return an expression of a matrix with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with entries equal to one
auto ones(Index rows, Index cols) -> decltype(Matrix::Ones(rows, cols));

/// Return an expression of a matrix with random entries
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with random entries
auto random(Index rows, Index cols) -> decltype(Matrix::Random(rows, cols));

/// Return an expression of an identity matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of an identity matrix
auto identity(Index rows, Index cols) -> decltype(Matrix::Identity(rows, cols));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename Derived>
auto rows(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename Derived>
auto rows(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num));

/// Return a view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived, typename Indices>
auto rows(Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> Eigen::MatrixRowsView<Derived, Indices>;

/// Return a const view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived, typename Indices>
auto rows(const Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> Eigen::MatrixRowsViewConst<Derived, Indices>;

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename Derived>
auto cols(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num));

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename Derived>
auto cols(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num));

/// Return a view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto cols(Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> Eigen::MatrixColsView<Derived, Indices>;

/// Return a const view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto cols(const Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> Eigen::MatrixColsViewConst<Derived, Indices>;

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto segment(Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto segment(const Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto block(Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto block(const Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto submatrix(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> Eigen::MatrixSubView<Derived, Indices>;

/// Return a const view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto submatrix(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> Eigen::MatrixSubViewConst<Derived, Indices>;

/// Return the transpose of the matrix
template<typename Derived>
auto tr(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose());

/// Return the transpose of the matrix
template<typename Derived>
auto tr(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose());

/// Return the component-wise inverse of the matrix
template<typename Derived>
auto inv(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseInverse());

/// Return a diagonal matrix representation of a vector
template<typename Derived>
auto diag(const Eigen::MatrixBase<Derived>& vec) -> decltype(vec.asDiagonal());

/// Return a vector representation of the diagonal of a matrix
template<typename Derived>
auto diagonal(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal());

/// Return a vector representation of the diagonal of a matrix
template<typename Derived>
auto diagonal(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal());

/// Return the Lp norm of a matrix
template<int p, typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the L2 norm of a matrix
template<typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the L-inf norm of a matrix
template<typename Derived>
auto norminf(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the sum of the components of a matrix
template<typename Derived>
auto sum(const Eigen::DenseBase<Derived>& mat) -> typename Derived::Scalar;

/// Return the dot product of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto dot(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.dot(rhs));

/// Return the minimum component of a matrix
template<typename Derived>
auto min(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.minCoeff());

/// Return the component-wise minimum of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto min(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMin(rhs));

/// Return the maximum component of a matrix
template<typename Derived>
auto max(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.maxCoeff());

/// Return the component-wise maximum of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto max(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMax(rhs));

/// Return the component-wise multiplication of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto operator%(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseProduct(rhs));

/// Return the component-wise division of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto operator/(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseQuotient(rhs));

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator/(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype(scalar*mat.cwiseInverse());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator+(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar + mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator+(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((scalar + mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator-(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar - mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator-(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((mat.array() - scalar).matrix());

/// Return the component-wise absolute entries of a matrix
template<typename Derived>
auto abs(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseAbs());

/// Return the component-wise square root of a matrix
template<typename Derived>
auto sqrt(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseSqrt());

/// Return the component-wise exponential of a matrix
template<typename Derived>
auto pow(const Eigen::MatrixBase<Derived>& mat, double power) -> decltype(mat.array().pow(power));

/// Return the component-wise natural exponential of a matrix
template<typename Derived>
auto exp(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().exp().matrix());

/// Return the component-wise natural log of a matrix
template<typename Derived>
auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log().matrix());

/// Return the component-wise log10 of a matrix
template<typename Derived>
auto log10(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log10().matrix());

} // namespace Reaktoro

#include "Matrix.hxx"
