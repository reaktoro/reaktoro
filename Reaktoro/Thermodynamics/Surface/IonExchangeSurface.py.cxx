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
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>

using namespace Reaktoro;

void exportIonExchangeSurface(py::module& m)
{
    py::class_<IonExchangeSurfaceState>(m, "IonExchangeSurfaceState")
    .def(py::init<>())
    .def_readwrite("meq", &IonExchangeSurfaceState::meq)
    .def_readwrite("lng", &IonExchangeSurfaceState::lng)
    .def_readwrite("n", &IonExchangeSurfaceState::n)
    .def_readwrite("meq", &IonExchangeSurfaceState::meq)
    ;

    py::class_<IonExchangeSurface>(m, "IonExchangeSurface")
    .def(py::init<const SpeciesList&>())
    .def("clone", &IonExchangeSurface::clone)
    .def("species", py::overload_cast<Index>(&IonExchangeSurface::species, py::const_))
    .def("species", py::overload_cast<>(&IonExchangeSurface::species, py::const_))
    .def("ze", &IonExchangeSurface::ze)
    .def("setLogarithmsOfActivities", &IonExchangeSurface::setLogarithmsOfActivities)
    .def("state", &IonExchangeSurface::state)
    ;
}
