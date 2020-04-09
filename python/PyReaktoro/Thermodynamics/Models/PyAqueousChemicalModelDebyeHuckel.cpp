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

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelDebyeHuckel.hpp>

namespace Reaktoro {

void exportAqueousChemicalModelDebyeHuckel(py::module& m)
{
    auto aiondefault1 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::aiondefault);
    auto aiondefault2 = static_cast<double(DebyeHuckelParams::*)() const>(&DebyeHuckelParams::aiondefault);

    auto aion1 = static_cast<void(DebyeHuckelParams::*)(std::string,double)>(&DebyeHuckelParams::aion);
    auto aion2 = static_cast<void(DebyeHuckelParams::*)(const std::map<std::string, double>&)>(&DebyeHuckelParams::aion);
    auto aion3 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::aion);
    auto aion4 = static_cast<double(DebyeHuckelParams::*)(std::string) const>(&DebyeHuckelParams::aion);

    auto biondefault1 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::biondefault);
    auto biondefault2 = static_cast<double(DebyeHuckelParams::*)() const>(&DebyeHuckelParams::biondefault);

    auto bion1 = static_cast<void(DebyeHuckelParams::*)(std::string,double)>(&DebyeHuckelParams::bion);
    auto bion2 = static_cast<void(DebyeHuckelParams::*)(const std::map<std::string, double>&)>(&DebyeHuckelParams::bion);
    auto bion3 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::bion);
    auto bion4 = static_cast<double(DebyeHuckelParams::*)(std::string) const>(&DebyeHuckelParams::bion);

    auto bneutraldefault1 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::bneutraldefault);
    auto bneutraldefault2 = static_cast<double(DebyeHuckelParams::*)() const>(&DebyeHuckelParams::bneutraldefault);

    auto bneutral1 = static_cast<void(DebyeHuckelParams::*)(std::string,double)>(&DebyeHuckelParams::bneutral);
    auto bneutral2 = static_cast<void(DebyeHuckelParams::*)(const std::map<std::string, double>&)>(&DebyeHuckelParams::bneutral);
    auto bneutral3 = static_cast<void(DebyeHuckelParams::*)(double)>(&DebyeHuckelParams::bneutral);
    auto bneutral4 = static_cast<double(DebyeHuckelParams::*)(std::string) const>(&DebyeHuckelParams::bneutral);

    py::class_<DebyeHuckelParams>(m, "DebyeHuckelParams")
        .def(py::init<>())
        .def("aiondefault", aiondefault1)
        .def("aiondefault", aiondefault2)
        .def("aion", aion1)
        .def("aion", aion2)
        .def("aion", aion3)
        .def("aion", aion4)
        .def("biondefault", biondefault1)
        .def("biondefault", biondefault2)
        .def("bion", bion1)
        .def("bion", bion2)
        .def("bion", bion3)
        .def("bion", bion4)
        .def("bneutraldefault", bneutraldefault1)
        .def("bneutraldefault", bneutraldefault2)
        .def("bneutral", bneutral1)
        .def("bneutral", bneutral2)
        .def("bneutral", bneutral3)
        .def("bneutral", bneutral4)
        .def("setLimitingLaw", &DebyeHuckelParams::setLimitingLaw)
        .def("setKielland1937", &DebyeHuckelParams::setKielland1937)
        .def("setWATEQ4F", &DebyeHuckelParams::setWATEQ4F)
        .def("setPHREEQC", &DebyeHuckelParams::setPHREEQC)
        ;
}

} // namespace Reaktoro
