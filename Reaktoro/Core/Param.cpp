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

#include "Param.hpp"

namespace Reaktoro {
namespace {

/// Issue a warning if @p val is out-of-bounds
auto warningIfOutOfBounds(const Param& param, const real& val) -> void
{
    const auto lb = param.lowerbound();
    const auto ub = param.upperbound();
    const auto id = param.id();
    warningif(val <= lb || val >= ub,
        "Setting parameter with value ", val, " and id `", id, "` violates "
        "either its lower bound (", lb, ") or its upper bound (", ub, ").")
}

} // namespace

/// The type containing the data members of the parameter.
struct Param::Impl
{
    /// The parameter value.
    real value = 0.0;

    /// The parameter id.
    String id;

    /// The parameter lower bound value (default: -inf).
    double lowerbound = -std::numeric_limits<double>::infinity();

    /// The parameter upper bound value (default: +inf).
    double upperbound = +std::numeric_limits<double>::infinity();

    /// The boolean flag that indicates if this parameter is constant.
    double isconst = false;
};

Param::Param()
: pimpl(new Impl())
{}

Param::Param(const real& val)
: pimpl(new Impl({ val }))
{}

Param::Param(double val)
: pimpl(new Impl({ val }))
{}

Param::Param(const String& id, const real& val)
: pimpl(new Impl({ val, id }))
{}

auto Param::clone() const -> Param
{
    Param param;
    *param.pimpl = *pimpl;
    return param;
}

auto Param::assign(const real& val) -> Param&
{
    warningIfOutOfBounds(*this, val);
    pimpl->value = val;
    return *this;
}

auto Param::value(const real& val) -> Param&
{
    return assign(val);;
}

auto Param::value() const -> const real&
{
    return pimpl->value;
}

auto Param::value() -> real&
{
    return pimpl->value;
}

auto Param::id(String id) -> Param&
{
    pimpl->id = id;
    return *this;
}

auto Param::id() const -> const String&
{
    return pimpl->id;
}

auto Param::lowerbound(double val) -> Param&
{
    pimpl->lowerbound = val;
    return *this;
}

auto Param::lowerbound() const -> double
{
    return pimpl->lowerbound;
}

auto Param::upperbound(double val) -> Param&
{
    pimpl->upperbound = val;
    return *this;
}

auto Param::upperbound() const -> double
{
    return pimpl->upperbound;
}

auto Param::isconst(bool val) -> Param&
{
    pimpl->isconst = val;
    return *this;
}

auto Param::isconst() const -> bool
{
    return pimpl->isconst;
}

auto Param::operator=(double val) -> Param&
{
    return value(val);
}

auto Param::operator=(const real& val) -> Param&
{
    return value(val);
}

Param::operator const real&() const
{
    return pimpl->value;
}

Param::operator real&()
{
    return pimpl->value;
}

Param::operator double() const
{
    return pimpl->value;
}

auto Param::Constant(const real& val) -> Param
{
    Param param(val);
    param.isconst(true);
    return param;
}

} // namespace Reaktoro
