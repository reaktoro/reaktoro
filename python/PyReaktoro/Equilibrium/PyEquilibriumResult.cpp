// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

void exportEquilibriumResult(py::module& m)
{
    py::class_<EquilibriumTiming>(m, "EquilibriumTiming")
        .def_readwrite("solve", &EquilibriumTiming::solve)
        .def_readwrite("standard_thermodynamic_properties", &EquilibriumTiming::standard_thermodynamic_properties)
        .def_readwrite("chemical_properties", &EquilibriumTiming::chemical_properties)
        ;

    py::class_<EquilibriumResult>(m, "EquilibriumResult")
        .def_readwrite("optimum", &EquilibriumResult::optimum)
        .def_readwrite("timing", &EquilibriumResult::timing)
        ;
}

} // namespace Reaktoro
