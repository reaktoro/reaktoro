// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// // pybind11 includes
// #include <pybind11/pybind11.h>
// #include <pybind11/functional.h>
// #include <pybind11/stl.h>
// namespace py = pybind11;

// // Reaktoro includes
// #include <Reaktoro/Core/ChemicalProps.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>
// #include <Reaktoro/Core/Species.hpp>
// #include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
// using namespace Reaktoro;

// void exportEquilibriumConstraints(py::module& m)
// {
//     const auto return_internal_ref = py::return_value_policy::reference_internal;

//     py::class_<EquilibriumConstraints>(m, "EquilibriumConstraints")
//         .def(py::init<const ChemicalSystem&>())
//         .def("control", &EquilibriumConstraints::control)
//         .def("until", &EquilibriumConstraints::until)
//         .def("preserve", &EquilibriumConstraints::preserve)
//         .def("fix", &EquilibriumConstraints::fix)
//         .def("prevent", &EquilibriumConstraints::prevent)
//         .def("system", &EquilibriumConstraints::system, return_internal_ref)
//         .def("data", &EquilibriumConstraints::data, return_internal_ref)
//         ;

//     py::class_<EquilibriumConstraints::Control>(m, "_EquilibriumConstraintsControl")
//         .def("temperature", &EquilibriumConstraints::Control::temperature, return_internal_ref)
//         .def("pressure", &EquilibriumConstraints::Control::pressure, return_internal_ref)
//         .def("titrationOf", &EquilibriumConstraints::Control::titrationOf, return_internal_ref)
//         .def("titrationOfEither", &EquilibriumConstraints::Control::titrationOfEither, return_internal_ref)
//         ;

//     py::class_<EquilibriumConstraints::Until>(m, "_EquilibriumConstraintsUntil")
//         .def("volume", &EquilibriumConstraints::Until::volume, return_internal_ref)
//         .def("internalEnergy", &EquilibriumConstraints::Until::internalEnergy, return_internal_ref)
//         .def("enthalpy", &EquilibriumConstraints::Until::enthalpy, return_internal_ref)
//         .def("gibbsEnergy", &EquilibriumConstraints::Until::gibbsEnergy, return_internal_ref)
//         .def("helmholtzEnergy", &EquilibriumConstraints::Until::helmholtzEnergy, return_internal_ref)
//         .def("entropy", &EquilibriumConstraints::Until::entropy, return_internal_ref)
//         .def("custom", &EquilibriumConstraints::Until::custom, return_internal_ref)
//         ;

//     py::class_<EquilibriumConstraints::Preserve>(m, "_EquilibriumConstraintsPreserve")
//         .def("volume", &EquilibriumConstraints::Preserve::volume, return_internal_ref)
//         .def("internalEnergy", &EquilibriumConstraints::Preserve::internalEnergy, return_internal_ref)
//         .def("enthalpy", &EquilibriumConstraints::Preserve::enthalpy, return_internal_ref)
//         .def("gibbsEnergy", &EquilibriumConstraints::Preserve::gibbsEnergy, return_internal_ref)
//         .def("helmholtzEnergy", &EquilibriumConstraints::Preserve::helmholtzEnergy, return_internal_ref)
//         .def("entropy", &EquilibriumConstraints::Preserve::entropy, return_internal_ref)
//         .def("custom", &EquilibriumConstraints::Preserve::custom, return_internal_ref)
//         ;

//     py::class_<EquilibriumConstraints::Fix>(m, "_EquilibriumConstraintsFix")
//         .def("chemicalPotential", py::overload_cast<const ChemicalFormula&, const Fn<real(real,real)>&>(&EquilibriumConstraints::Fix::chemicalPotential))
//         .def("chemicalPotential", py::overload_cast<String, real, String>(&EquilibriumConstraints::Fix::chemicalPotential))
//         .def("lnActivity", py::overload_cast<const Species&, real>(&EquilibriumConstraints::Fix::lnActivity))
//         .def("lnActivity", py::overload_cast<String, real>(&EquilibriumConstraints::Fix::lnActivity))
//         .def("lgActivity", &EquilibriumConstraints::Fix::lgActivity, return_internal_ref)
//         .def("activity", &EquilibriumConstraints::Fix::activity, return_internal_ref)
//         .def("fugacity", &EquilibriumConstraints::Fix::fugacity, return_internal_ref)
//         .def("pH", &EquilibriumConstraints::Fix::pH, return_internal_ref)
//         .def("pMg", &EquilibriumConstraints::Fix::pMg, return_internal_ref)
//         .def("pe", &EquilibriumConstraints::Fix::pe, return_internal_ref)
//         .def("Eh", &EquilibriumConstraints::Fix::Eh, return_internal_ref)
//         ;

//     py::class_<EquilibriumConstraints::Prevent>(m, "_EquilibriumConstraintsPrevent")
//         .def("fromReacting", py::overload_cast<Index>(&EquilibriumConstraints::Prevent::fromReacting), return_internal_ref)
//         .def("fromReacting", py::overload_cast<Pairs<Index, double>>(&EquilibriumConstraints::Prevent::fromReacting), return_internal_ref)
//         .def("fromReacting", py::overload_cast<String>(&EquilibriumConstraints::Prevent::fromReacting), return_internal_ref)
//         .def("fromIncreasing", py::overload_cast<Index>(&EquilibriumConstraints::Prevent::fromIncreasing), return_internal_ref)
//         .def("fromIncreasing", py::overload_cast<String>(&EquilibriumConstraints::Prevent::fromIncreasing), return_internal_ref)
//         .def("fromDecreasing", py::overload_cast<Index>(&EquilibriumConstraints::Prevent::fromDecreasing), return_internal_ref)
//         .def("fromDecreasing", py::overload_cast<String>(&EquilibriumConstraints::Prevent::fromDecreasing), return_internal_ref)
//         ;

//     py::class_<EquilibriumEquationArgs>(m, "EquilibriumEquationArgs")
//         .def_property_readonly("props", [](const EquilibriumEquationArgs& self) { return self.props; })
//         .def_readonly("q", &EquilibriumEquationArgs::q)
//         .def_property_readonly("titrants", [](const EquilibriumEquationArgs& self) { return self.titrants; })
//         .def("titrantAmount", &EquilibriumEquationArgs::titrantAmount)
//         ;
// }
