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
        .def_readwrite("solve", &SmartEquilibriumTiming::solve, "The time spent in seconds for solving the chemical equilibrium problem.")
        .def_readwrite("learn", &SmartEquilibriumTiming::learn, "The time spent in seconds for learning a new chemical equilibrium calculation.")
        .def_readwrite("learn_gibbs_energy_minimization", &SmartEquilibriumTiming::learn_gibbs_energy_minimization, "The time spent in seconds for a conventional Gibbs energy minimization calculation during learning operation.")
        .def_readwrite("learn_chemical_properties", &SmartEquilibriumTiming::learn_chemical_properties, "The time spent in seconds for computing the chemical properties during learning operation.")
        .def_readwrite("learn_sensitivity_matrix", &SmartEquilibriumTiming::learn_sensitivity_matrix, "The time spent in seconds for computing the sensitivity matrix during learning operation.")
        .def_readwrite("learn_error_control_matrices", &SmartEquilibriumTiming::learn_error_control_matrices, "The time spent in seconds for computing the error control matrices during learning operation.")
        .def_readwrite("learn_storage", &SmartEquilibriumTiming::learn_storage, "The time spent in seconds for storing the computed chemical state into the tree of knowledge.")
        .def_readwrite("estimate", &SmartEquilibriumTiming::estimate, "The time spent in seconds for the smart chemical equilibrium state estimation.")
        .def_readwrite("estimate_search", &SmartEquilibriumTiming::estimate_search, "The time spent in seconds for the search operation during a smart estimation.")
        .def_readwrite("estimate_error_control", &SmartEquilibriumTiming::estimate_error_control, "The time spent in seconds during on error control while searching during a smart estimation.")
        .def_readwrite("estimate_taylor", &SmartEquilibriumTiming::estimate_taylor, "The time spent in seconds for the matrix-vector multiplication during a smart estimation.")
        .def_readwrite("estimate_database_priority_update", &SmartEquilibriumTiming::estimate_database_priority_update, "The time spent in seconds for updating the priority related info of the clusters in the database.")
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
