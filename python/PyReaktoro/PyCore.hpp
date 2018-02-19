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
#include <PyReaktoro/Core/PyChemicalOutput.hpp>
#include <PyReaktoro/Core/PyChemicalPlot.hpp>
#include <PyReaktoro/Core/PyChemicalProperties.hpp>
#include <PyReaktoro/Core/PyChemicalProperty.hpp>
#include <PyReaktoro/Core/PyChemicalQuantity.hpp>
#include <PyReaktoro/Core/PyChemicalState.hpp>
#include <PyReaktoro/Core/PyChemicalSystem.hpp>
#include <PyReaktoro/Core/PyConnectivity.hpp>
#include <PyReaktoro/Core/PyElement.hpp>
#include <PyReaktoro/Core/PyPartition.hpp>
#include <PyReaktoro/Core/PyPhase.hpp>
#include <PyReaktoro/Core/PyReaction.hpp>
#include <PyReaktoro/Core/PyReactionSystem.hpp>
#include <PyReaktoro/Core/PySpecies.hpp>
#include <PyReaktoro/Core/PyThermoProperties.hpp>
#include <PyReaktoro/Core/PyUtils.hpp>

namespace Reaktoro {

inline auto export_Core() -> void
{
    export_ChemicalOutput();
    export_ChemicalPlot();
    export_ChemicalProperties();
    export_ChemicalProperty();
    export_ChemicalQuantity();
    export_ChemicalState();
    export_ChemicalSystem();
    export_Connectivity();
    export_Element();
    export_Partition();
    export_Phase();
    export_Reaction();
    export_ReactionSystem();
    export_Species();
    export_ThermoProperties();
    export_Utils();
}

} // namespace Reaktoro
