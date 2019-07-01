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
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {

void exportSmartEquilibriumResult(py::module& m)
{
    py::class_<SmartEquilibriumTiming>(m, "SmartEquilibriumTiming")
        .def_readwrite("solve", &SmartEquilibriumTiming::solve)
        .def_readwrite("learning", &SmartEquilibriumTiming::learning)
        .def_readwrite("learning_gibbs_energy_minimization", &SmartEquilibriumTiming::learning_gibbs_energy_minimization)
        .def_readwrite("learning_chemical_properties", &SmartEquilibriumTiming::learning_chemical_properties)
        .def_readwrite("learning_sensitivity_matrix", &SmartEquilibriumTiming::learning_sensitivity_matrix)
        .def_readwrite("learning_storage", &SmartEquilibriumTiming::learning_storage)
        .def_readwrite("estimate", &SmartEquilibriumTiming::estimate)
        .def_readwrite("estimate_search", &SmartEquilibriumTiming::estimate_search)
        .def_readwrite("estimate_mat_vec_mul", &SmartEquilibriumTiming::estimate_mat_vec_mul)
        .def_readwrite("estimate_acceptance", &SmartEquilibriumTiming::estimate_acceptance)
        .def(py::self += py::self)
        ;

    py::class_<SmartEquilibriumResultDuringEstimate>(m, "SmartEquilibriumResultDuringEstimate")
        .def_readwrite("accepted", &SmartEquilibriumResultDuringEstimate::accepted)
        .def_readwrite("failed_with_species", &SmartEquilibriumResultDuringEstimate::failed_with_species)
        .def_readwrite("failed_with_amount", &SmartEquilibriumResultDuringEstimate::failed_with_amount)
        .def_readwrite("failed_with_chemical_potential", &SmartEquilibriumResultDuringEstimate::failed_with_chemical_potential)
        ;

    py::class_<SmartEquilibriumResultDuringLearning>(m, "SmartEquilibriumResultDuringLearning")
        .def_readwrite("gibbs_energy_minimization", &SmartEquilibriumResultDuringLearning::gibbs_energy_minimization)
        ;

    py::class_<SmartEquilibriumResult>(m, "SmartEquilibriumResult")
        .def_readwrite("estimate", &SmartEquilibriumResult::estimate)
        .def_readwrite("learning", &SmartEquilibriumResult::learning)
        .def_readwrite("timing", &SmartEquilibriumResult::timing)
        ;
}

} // namespace Reaktoro
