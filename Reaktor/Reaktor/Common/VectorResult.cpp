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

#include "VectorResult.hpp"

// Reaktor includes
#include <Reaktor/Common/ScalarResult.hpp>

namespace Reaktor {

VectorResult::VectorResult()
{}

VectorResult::VectorResult(unsigned dim)
: VectorResult(dim, dim)
{}

VectorResult::VectorResult(unsigned domain, unsigned range)
: val(Vector(range)), grad(Matrix(range, domain))
{}

VectorResult::VectorResult(const Vector& val, const Matrix& grad)
: val(val), grad(grad)
{}

auto VectorResult::row(const Index& i) -> VectorResultRow
{
	return VectorResultRow(val.row(i), grad.row(i));
}

auto VectorResult::row(const Index& i) const -> const VectorResultRow
{
	return VectorResultRow(val.row(i), grad.row(i));
}

auto VectorResultRow::operator=(const ScalarResult& res) -> VectorResultRow&
{
	val = res.val;
	grad = res.grad;
	return *this;
}

} // namespace Reaktor
