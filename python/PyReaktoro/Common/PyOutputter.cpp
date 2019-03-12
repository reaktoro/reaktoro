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
#include <Reaktoro/Common/Outputter.hpp>

namespace Reaktoro {

void exportOutputter(py::module& m)
{
    py::class_<OutputterOptions>(m, "OutputterOptions")
        .def_readwrite("active", &OutputterOptions::active)
        .def_readwrite("fixed", &OutputterOptions::fixed)
        .def_readwrite("scientific", &OutputterOptions::scientific)
        .def_readwrite("precision", &OutputterOptions::precision)
        .def_readwrite("width", &OutputterOptions::width)
        .def_readwrite("separator", &OutputterOptions::separator)
        ;
}

} // namespace Reaktoro

