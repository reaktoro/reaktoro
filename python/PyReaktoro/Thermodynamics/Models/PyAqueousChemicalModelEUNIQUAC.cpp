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
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelEUNIQUAC.hpp>

namespace Reaktoro {

void exportAqueousChemicalModelEUNIQUAC(py::module& m)
{
    auto ri1 = static_cast<void(EUNIQUACParams::*)(const std::string&,double)>(&EUNIQUACParams::ri);
    auto ri2 = static_cast<void(EUNIQUACParams::*)(const std::map<std::string, double>&)>(&EUNIQUACParams::ri);
    auto ri3 = static_cast<double(EUNIQUACParams::*)(const std::string&) const>(&EUNIQUACParams::ri);

    py::class_<EUNIQUACParams>(m, "EUNIQUACParams")
        .def(py::init<>())
        .def("ri", ri1)
        .def("ri", ri2)
        .def("ri", ri3)
        .def("setDTUvalues", &EUNIQUACParams::setDTUvalues)
        ;
}

} // namespace Reaktoro
