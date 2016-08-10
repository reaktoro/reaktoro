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
#include <PyReaktoro/Thermodynamics/Phases/PyAqueousPhase.hpp>
#include <PyReaktoro/Thermodynamics/Phases/PyGaseousPhase.hpp>
#include <PyReaktoro/Thermodynamics/Phases/PyMineralPhase.hpp>
#include <PyReaktoro/Thermodynamics/PyChemicalEditor.hpp>
#include <PyReaktoro/Thermodynamics/PyDatabase.hpp>
#include <PyReaktoro/Thermodynamics/PyThermo.hpp>
#include <PyReaktoro/Thermodynamics/Species/PyAqueousSpecies.hpp>
#include <PyReaktoro/Thermodynamics/Species/PyGaseousSpecies.hpp>
#include <PyReaktoro/Thermodynamics/Species/PyMineralSpecies.hpp>

namespace Reaktoro {

inline auto export_Thermodynamics() -> void
{
    export_Database();
    export_ChemicalEditor();
    export_Thermo();
    export_AqueousPhase();
    export_GaseousPhase();
    export_MineralPhase();
    export_AqueousSpecies();
    export_GaseousSpecies();
    export_MineralSpecies();
}

} // namespace Reaktoro
