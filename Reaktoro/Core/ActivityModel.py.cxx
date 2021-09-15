// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Model.py.hxx>
using namespace Reaktoro;

void exportActivityModel(py::module& m)
{
    py::class_<ActivityArgs>(m, "ActivityArgs")
        .def(py::init<const real&, const real&, ArrayXrConstRef>())
        .def_property_readonly("T", [](const ActivityArgs& self) { return self.T; })
        .def_property_readonly("P", [](const ActivityArgs& self) { return self.P; })
        .def_property_readonly("x", [](const ActivityArgs& self) { return self.x; })
        ;

    auto cls = exportModel<ActivityProps, ActivityArgs>(m, "ActivityModel");

    cls.def("__call__", [](const ActivityModel& self, const real& T, const real& P, ArrayXrConstRef x) { return self({T, P, x}); });
    cls.def("__call__", [](const ActivityModel& self, ActivityPropsRef props, const real& T, const real& P, ArrayXrConstRef x) { self(props, {T, P, x}); });

    auto chain4py = [](py::args args)
    {
        Vec<ActivityModelGenerator> models;
        models.reserve(args.size());
        for(auto arg : args)
            models.push_back(arg.cast<ActivityModelGenerator>());
        return chain(models);
    };

    m.def("chain", chain4py);
}
