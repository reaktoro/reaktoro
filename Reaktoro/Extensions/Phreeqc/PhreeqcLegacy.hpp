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

#pragma once

// Phreeqc includes
#define Phreeqc PHREEQC
#define protected public

#include <phreeqc/Phreeqc.h>
#include <phreeqc/GasPhase.h>

using PhreeqcElement = element;
using PhreeqcSpecies = species;
using PhreeqcPhase = phase;

#undef Phreeqc
#undef protected
#undef pi

//==================================================
// WARNING WARNING WARNING WARNING WARNING WARNING
//==================================================
// This header file must not be included by another
// header file that will be exposed to users.
// If so, this propagates the need for phreeqc
// header files to be available in the user system.
//==================================================
