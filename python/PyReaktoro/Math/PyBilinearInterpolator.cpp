// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {

void exportBilinearInterpolator(py::module& m)
{
    py::class_<BilinearInterpolator>(m, "BilinearInterpolator")
        .def(py::init<>())
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const std::vector<double>&>())
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const std::function<double(double, double)>&>())
        .def("setCoordinatesX", &BilinearInterpolator::setCoordinatesX)
        .def("setCoordinatesY", &BilinearInterpolator::setCoordinatesY)
        .def("setData", &BilinearInterpolator::setData)
        .def("xCoordinates", &BilinearInterpolator::xCoordinates, py::return_value_policy::reference_internal)
        .def("yCoordinates", &BilinearInterpolator::yCoordinates, py::return_value_policy::reference_internal)
        .def("data", &BilinearInterpolator::data, py::return_value_policy::reference_internal)
        .def("empty", &BilinearInterpolator::empty)
        .def("__call__", &BilinearInterpolator::operator())
        ;
}

} // namespace Reaktoro
