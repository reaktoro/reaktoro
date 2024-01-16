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
#include <Reaktoro/Utils/Material.hpp>
using namespace Reaktoro;

void exportMaterial(py::module& m)
{
    py::class_<Material>(m, "Material")
        .def(py::init<ChemicalSystem const&>())
        .def("addSpeciesAmount", py::overload_cast<StringOrIndex const&, double>(&Material::addSpeciesAmount))
        .def("addSpeciesAmount", py::overload_cast<StringOrIndex const&, double, Chars>(&Material::addSpeciesAmount))
        .def("addSpeciesMass", &Material::addSpeciesMass)
        .def("addSubstanceAmount", py::overload_cast<ChemicalFormula const&, double>(&Material::addSubstanceAmount))
        .def("addSubstanceAmount", py::overload_cast<ChemicalFormula const&, double, Chars>(&Material::addSubstanceAmount))
        .def("addSubstanceMass", &Material::addSubstanceMass)
        .def("addMaterialAmount", py::overload_cast<Material const&, double>(&Material::addMaterialAmount))
        .def("addMaterialAmount", py::overload_cast<Material const&, double, Chars>(&Material::addMaterialAmount))
        .def("addMaterialMass", &Material::addMaterialMass)
        .def("add", py::overload_cast<String const&, double, Chars>(&Material::add))
        .def("add", py::overload_cast<Material const&, double, Chars>(&Material::add))
        .def("scaleAmount", &Material::scaleAmount)
        .def("scaleMass", &Material::scaleMass)
        .def("scale", &Material::scale)
        // .def("with", &Material::with)  // with is a keyword in Python
        .def("system", &Material::system)
        .def("substances", &Material::substances, return_internal_ref)
        .def("species", &Material::species, return_internal_ref)
        .def("componentAmounts", &Material::componentAmounts)
        .def("elementAmounts", &Material::elementAmounts)
        .def("charge", &Material::charge)
        .def("amount", &Material::amount)
        .def("mass", &Material::mass)
        .def("molarMass", &Material::molarMass)
        .def("equilibrate", py::overload_cast<>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<EquilibriumOptions const&>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<EquilibriumRestrictions const&>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<EquilibriumRestrictions const&, EquilibriumOptions const&>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<double, Chars, double, Chars>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<double, Chars, double, Chars, EquilibriumOptions const&>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<double, Chars, double, Chars, EquilibriumRestrictions const&>(&Material::equilibrate))
        .def("equilibrate", py::overload_cast<double, Chars, double, Chars, EquilibriumRestrictions const&, EquilibriumOptions const&>(&Material::equilibrate))
        .def("result", &Material::result, return_internal_ref)
        .def("initialState", &Material::initialState)
        .def("__call__", [](const Material& self, double value, Chars unit) { return self(value, unit); })
        .def("__repr__", [](const Material& self) { std::stringstream ss; ss << self; return ss.str(); })
        .def(py::self + py::self)
        ;
}
