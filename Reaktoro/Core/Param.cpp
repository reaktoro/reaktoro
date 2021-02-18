// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
    const auto name = param.name();
    warningif(val <= lb || val >= ub,
        "Setting parameter with value ", val, " and name `", name, "` violates "
        "either its lower bound (", lb, ") or its upper bound (", ub, ").")
}

} // namespace

/// The type containing the data members of the parameter.
struct Param::Impl
{
    /// The parameter value.
    real value = 0.0;

    /// The parameter name.
    String name;

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

auto Param::value(const real& val) -> Param&
{
    warningIfOutOfBounds(*this, val);
    pimpl->value = val;
    return *this;
}

auto Param::value() const -> const real&
{
    return pimpl->value;
}

auto Param::name(String name) -> Param&
{
    pimpl->name = name;
    return *this;
}

auto Param::name() const -> const String&
{
    return pimpl->name;
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

auto Param::operator=(const real& val) -> Param&
{
    value(val);
    return *this;
}

Param::operator const real&() const
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
