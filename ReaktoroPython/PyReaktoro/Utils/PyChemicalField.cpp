// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Utils/ChemicalField.hpp>

// PyReaktoro includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

/// Dummy comparison to satisfy requirements of Boost.Python std::vector wrapper
auto operator<(const ChemicalField& lhs, const ChemicalField& rhs) -> bool
{
    return lhs.size() < rhs.size();
}

/// Dummy comparison to satisfy requirements of Boost.Python std::vector wrapper
auto operator==(const ChemicalField& lhs, const ChemicalField& rhs) -> bool
{
    return lhs.size() == rhs.size();
}

auto export_ChemicalField() -> void
{
    py::class_<ChemicalField>("ChemicalField")
        .def(py::init<>())
        .def(py::init<const Partition&, Index>())
        .def("set", &ChemicalField::set)
        .def("size", &ChemicalField::size)
        .def_readwrite("val", &ChemicalField::val)
        .def_readwrite("T", &ChemicalField::T)
        .def_readwrite("P", &ChemicalField::P)
        .def_readwrite("be", &ChemicalField::be)
        .def_readwrite("nk", &ChemicalField::nk)
        ;

    export_std_vector<ChemicalField>("ChemicalFieldVector");
}

} // namespace Reaktoro
