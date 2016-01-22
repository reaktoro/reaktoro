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

// PyReaktoro includes
#include <PyReaktoro/Interfaces/PyGems.hpp>
#include <PyReaktoro/Interfaces/PyInterface.hpp>
#include <PyReaktoro/Interfaces/PyPhreeqc.hpp>

namespace Reaktoro {

inline auto export_Interfaces() -> void
{
    // Warning: export_Interface() should always come first, otherwise the following
    // RuntimeError will be raised in python: "RuntimeError: extension class wrapper
    // for base class has not been created yet".
    export_Interface();
    export_Gems();
    export_Phreeqc();
}

} // namespace Reaktoro
