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
//
//#include "PyChemicalField.hpp"
//
//// Boost includes
//#include <boost/python.hpp>
//namespace py = boost::python;
//
//// Reaktoro includes
//#include <Reaktoro/Core/Partition.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
//#include <Reaktoro/Util/ChemicalField.hpp>
//
//// PyReaktoro includes
//#include <PyReaktoro/Common/PyConverters.hpp>
//
//namespace Reaktoro {
//
///// Dummy comparison to satisfy requirements of Boost.Python std::vector wrapper
//auto operator<(const ChemicalField& lhs, const ChemicalField& rhs) -> bool
//{
//    return lhs.size() < rhs.size();
//}
//
///// Dummy comparison to satisfy requirements of Boost.Python std::vector wrapper
//auto operator==(const ChemicalField& lhs, const ChemicalField& rhs) -> bool
//{
//    return lhs.size() == rhs.size();
//}
//
//void exportChemicalField(py::module& m)
//{
//    auto val = static_cast<VectorRef(ChemicalField::*)()>(&ChemicalField::val);
//    auto ddT = static_cast<VectorRef(ChemicalField::*)()>(&ChemicalField::ddT);
//    auto ddP = static_cast<VectorRef(ChemicalField::*)()>(&ChemicalField::ddP);
//    auto ddbe = static_cast<std::vector<Vector>&(ChemicalField::*)()>(&ChemicalField::ddbe);
//    auto ddnk = static_cast<std::vector<Vector>&(ChemicalField::*)()>(&ChemicalField::ddnk);
//
//    py::class_<ChemicalField>(m, "ChemicalField")
//        .def(py::init<>())
//        .def(py::init<const Partition&, Index>())
//        .def("set", &ChemicalField::set)
//        .def("partition", &ChemicalField::partition, py::return_value_policy::reference_internal)
//        .def("size", &ChemicalField::size)
//        .def("val", val, py::return_value_policy::reference_internal)
//        .def("ddT", ddT, py::return_value_policy::reference_internal)
//        .def("ddP", ddP, py::return_value_policy::reference_internal)
//        .def("ddbe", ddbe, py::return_value_policy::reference_internal)
//        .def("ddnk", ddnk, py::return_value_policy::reference_internal)
//        .def(py::self_ns::str(py::self_ns::self))
//        ;
//
//    exportstd_vector<ChemicalField>("ChemicalFieldVector");
//}
//
//} // namespace Reaktoro
