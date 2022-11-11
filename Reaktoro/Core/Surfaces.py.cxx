// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Surface.hpp>
#include <Reaktoro/Core/Surfaces.hpp>
using namespace Reaktoro;

void exportSurfaces(py::module& m)
{
    py::class_<GeneralSurface>(m, "GeneralSurface")
        .def(py::init<>())
        .def(py::init<String const&>())
        .def("setName", &GeneralSurface::setName, "Set the unique name of the surface.")
        .def("setAreaModel", &GeneralSurface::setAreaModel, "Set the area model of the surface.")
        .def("set", &GeneralSurface::set, "Set the area model of the surface (equivalent to GeneralSurface::setAreaModel).")
        .def("name", &GeneralSurface::name, "Return the name of the surface.")
        .def("areaModel", &GeneralSurface::areaModel, "Return the area model of the surface.")
        .def("convert", &GeneralSurface::operator(), "Convert this GeneralSurface object into a Surface object.") // NOTE: Do not use __call__ here because pybind11 will gladly cast a Python GeneralSurface object to a std::function of any type without any runtime errors! When checking if an argument in a ChemicalSystem constructor is of type ReactionGenerator or SurfaceGenerator (both objects of class std::function), the Python GeneralSurface object will be sucessfully converted, which is not expected.
        ;

    auto add = [](Surfaces& self, py::object generator)
    {
        try { self.add(generator.cast<Surface const&>()); }
        catch(...) {
            try { self.add(generator.cast<GeneralSurface const&>()); }
            catch(...) {
                try { self.add(generator.cast<SurfaceGenerator>()); }
                catch(...) {
                    errorif(true, "Could not add reaction generator object to Surfaces object:\n", py::str(generator));
                }
            }
        }
    };

    auto createSurfaces = [](py::args reaction_generators)
    {
        Surfaces reactions;
        for(auto generator : reaction_generators)
        {
            try { reactions.add(generator.cast<Surface const&>()); }
            catch(...) {
                try { reactions.add(generator.cast<GeneralSurface const&>()); }
                catch(...) {
                    try { reactions.add(generator.cast<SurfaceGenerator>()); }
                    catch(...) {
                        errorif(true, "Could not create Surfaces with reaction generator object:\n", py::str(generator));
                    }
                }
            }
        }

        return reactions;
    };

    py::class_<Surfaces>(m, "Surfaces")
        .def(py::init<>())
        .def(py::init(createSurfaces))
        .def("add", add)
        .def("convert", &Surfaces::convert)
        ;
}
