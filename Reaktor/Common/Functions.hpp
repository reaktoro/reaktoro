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

#pragma once

// C++ includes
#include <functional>

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ThermoScalar;
class ThermoVector;

typedef std::function<double(double, double)> FunctionTP;

typedef std::function<ThermoScalar(double, double, const Vector&)> ScalarFunctionTPn;

typedef std::function<ThermoVector(double, double, const Vector&)> VectorFunctionTPn;

typedef std::function<ThermoScalar(double, double, const Vector&, const ThermoVector&)> ScalarFunctionTPna;

typedef std::function<ThermoVector(double, double, const Vector&, const ThermoVector&)> VectorFunctionTPna;

/// Return a constant FunctionTP instance
inline auto constantFunctionTP(double value) -> FunctionTP
{
	return [=](double, double) { return value; };
}

} // namespace Reaktor
