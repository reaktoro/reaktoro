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
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>

namespace Reaktoro {

void exportSmartEquilibriumOptions(py::module& m)
{
    py::enum_<SmartEquilibriumStrategy>(m, "SmartEquilibriumStrategy")
        .value("NearestNeighbour", SmartEquilibriumStrategy::NearestNeighbour)
        .value("PriorityQueue", SmartEquilibriumStrategy::PriorityQueue)
        .value("Clustering", SmartEquilibriumStrategy::Clustering)
        ;

    py::class_<SmartEquilibriumOptions>(m, "SmartEquilibriumOptions")
        .def(py::init<>())
        .def_readwrite("learning", &SmartEquilibriumOptions::learning)
        .def_readwrite("mole_fraction_cutoff", &SmartEquilibriumOptions::mole_fraction_cutoff)
        .def_readwrite("amount_fraction_cutoff", &SmartEquilibriumOptions::amount_fraction_cutoff)
        .def_readwrite("reltol", &SmartEquilibriumOptions::reltol)
        .def_readwrite("method", &SmartEquilibriumOptions::method)
        ;
}

} // namespace Reaktoro
