// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/ReactionRate.hpp>
using namespace Reaktoro;

void exportReactionRate(py::module& m)
{
    py::class_<ReactionRate>(m, "ReactionRate")
        .def(py::init<>(), "Construct a default ReactionRate object.")
        .def(py::init<double>(), "Construct a ReactionRate object with given rate value.")
        .def(py::init<real const&>(), "Construct a ReactionRate object with given rate value.")
        .def_static("enforce", &ReactionRate::enforce, "Return a ReactionRate object that represents the residual of an enforced equation `f(props) = 0` instead of a reaction rate.")
        .def("value", &ReactionRate::value, return_internal_ref, "Return the underlying real object in the ReactionRate object.")
        .def("onEquationMode", &ReactionRate::onEquationMode, return_internal_ref, "Return true if this ReactionRate object is in equation mode.")
        .def("assign", [](ReactionRate& self, real const& value) { return self = value; }, return_internal_ref, "Assign a real value to this ReactionRate object.")

        .def(py::self += double())
        .def(py::self -= double())
        .def(py::self *= double())
        .def(py::self /= double())

        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(py::self / double())

        .def(double() + py::self)
        .def(double() - py::self)
        .def(double() * py::self)
        .def(double() / py::self)

        .def(py::self += real())
        .def(py::self -= real())
        .def(py::self *= real())
        .def(py::self /= real())

        .def(py::self + real())
        .def(py::self - real())
        .def(py::self * real())
        .def(py::self / real())

        .def(real() + py::self)
        .def(real() - py::self)
        .def(real() * py::self)
        .def(real() / py::self)

        .def(+py::self)
        .def(-py::self)
        ;

    // In case the user forgets to define a reaction rate function returning a
    // ReactionRate object, make sure the returning autodiff::real object can be
    // implictly convertible to the expected ReactionRate object.
    py::implicitly_convertible<real, ReactionRate>();
}
