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
#include <Reaktoro/Equilibrium/SmartEquilibriumProfiling.hpp>

namespace Reaktoro {

void exportSmartEquilibriumProfiling(py::module& m)
{
    py::class_<SmartEquilibriumProfiling>(m, "SmartEquilibriumProfiling")
        .def_readwrite("time_solve", &SmartEquilibriumProfiling::time_solve)
        .def_readwrite("time_learning", &SmartEquilibriumProfiling::time_learning)
        .def_readwrite("time_learning_gibbs_energy_minimization", &SmartEquilibriumProfiling::time_learning_gibbs_energy_minimization)
        .def_readwrite("time_learning_chemical_properties", &SmartEquilibriumProfiling::time_learning_chemical_properties)
        .def_readwrite("time_learning_sensitivity_matrix", &SmartEquilibriumProfiling::time_learning_sensitivity_matrix)
        .def_readwrite("time_learning_storage", &SmartEquilibriumProfiling::time_learning_storage)
        .def_readwrite("time_estimate", &SmartEquilibriumProfiling::time_estimate)
        .def_readwrite("time_estimate_search", &SmartEquilibriumProfiling::time_estimate_search)
        .def_readwrite("time_estimate_mat_vec_mul", &SmartEquilibriumProfiling::time_estimate_mat_vec_mul)
        .def_readwrite("time_estimate_acceptance", &SmartEquilibriumProfiling::time_estimate_acceptance)
        ;
}

} // namespace Reaktoro
