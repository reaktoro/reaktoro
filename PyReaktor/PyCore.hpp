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

// PyReaktor includes
#include <PyReaktor/Core/PyChemicalState.hpp>
#include <PyReaktor/Core/PyChemicalSystem.hpp>
#include <PyReaktor/Core/PyElement.hpp>
#include <PyReaktor/Core/PyPartition.hpp>
#include <PyReaktor/Core/PyPhase.hpp>
#include <PyReaktor/Core/PyReaction.hpp>
#include <PyReaktor/Core/PySpecies.hpp>
#include <PyReaktor/Core/PyUtils.hpp>

namespace Reaktor {

inline auto export_Core() -> void
{
    export_Element();
    export_Species();
    export_Phase();
    export_ChemicalSystem();
    export_ChemicalState();
    export_Partition();
    export_Reaction();
    export_Utils();
}

} // namespace Reaktor
