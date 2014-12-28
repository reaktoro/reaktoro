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

#include "PyOutputter.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {

auto export_Outputter() -> void
{
    py::class_<OutputOptions>("OutputOptions")
        .def_readwrite("active", &OutputOptions::active)
        .def_readwrite("fixed", &OutputOptions::fixed)
        .def_readwrite("scientific", &OutputOptions::scientific)
        .def_readwrite("precision", &OutputOptions::precision)
        .def_readwrite("width", &OutputOptions::width)
        .def_readwrite("separator", &OutputOptions::separator)
        ;
}

} // namespace Reaktor

