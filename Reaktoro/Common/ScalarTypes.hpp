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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

class Temperature;
class Pressure;

template<typename V>
class ThermoScalarBase;

template<typename V, typename T, typename P>
class ThermoVectorBase;

template<typename V, typename N>
class ChemicalScalarBase;

template<typename V, typename T, typename P, typename N>
class ChemicalVectorBase;

using ThermoScalar = ThermoScalarBase<double>;
using ThermoVector = ThermoVectorBase<Vector,Vector,Vector>;

using ChemicalScalar = ChemicalScalarBase<double,Vector>;
using ChemicalVector = ChemicalVectorBase<Vector,Vector,Vector,Matrix>;

using ThermoScalarFunction = std::function<ThermoScalar(Temperature, Pressure)>;
using ThermoVectorFunction = std::function<ThermoVector(Temperature, Pressure)>;

using ChemicalScalarFunction = std::function<ChemicalScalar(Temperature, Pressure, const Vector&)>;
using ChemicalVectorFunction = std::function<ChemicalVector(Temperature, Pressure, const Vector&)>;

} // namespace Reaktoro
