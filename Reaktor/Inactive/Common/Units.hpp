/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Units++ includes
#include <Reaktor/External/Units/Units.hpp>

namespace units {

using Density               = Constant<decltype(kg()/m3())>;
using MolarMass             = Constant<decltype(g()/mol())>;
using MolarVolume           = Constant<decltype(m3()/mol())>;
using SpecificSurfaceArea   = Constant<decltype(m2()/g())>;
using VolumetricSurfaceArea = Constant<decltype(m2()/m3())>;

} /* namespace units */

namespace Reaktor {

} // namespace Reaktor
