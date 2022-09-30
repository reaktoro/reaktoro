// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
using namespace Reaktoro;

void exportSmartEquilibriumResult(py::module& m)
{
    py::class_<SmartEquilibriumTiming>(m, "SmartEquilibriumTiming")
        .def_readwrite("solve", &SmartEquilibriumTiming::solve, "The time spent for solving the chemical equilibrium problem (in seconds).")
        .def_readwrite("learning", &SmartEquilibriumTiming::learning, "The time spent for learning a new chemical equilibrium state (in seconds).")
        .def_readwrite("learning_solve", &SmartEquilibriumTiming::learning_solve, "The time spent for a conventional iterative chemical equilibrium calculation during the learning operation (in seconds).")
        .def_readwrite("learning_chemical_properties", &SmartEquilibriumTiming::learning_chemical_properties, "The time spent for computing the chemical properties during the learning operation (in seconds).")
        .def_readwrite("learning_sensitivity_matrix", &SmartEquilibriumTiming::learning_sensitivity_matrix, "The time spent for computing the sensitivity matrix during the learning operation (in seconds).")
        .def_readwrite("learning_error_control_matrices", &SmartEquilibriumTiming::learning_error_control_matrices, "The time spent for computing the error control matrices during the learning operation (in seconds).")
        .def_readwrite("learning_storage", &SmartEquilibriumTiming::learning_storage, "The time spent for storing the computed chemical state into the tree of knowledge (in seconds).")
        .def_readwrite("prediction", &SmartEquilibriumTiming::prediction, "The time spent for the smart chemical equilibrium state prediction (in seconds).")
        .def_readwrite("prediction_search", &SmartEquilibriumTiming::prediction_search, "The time spent for the search operation during a smart prediction (in seconds).")
        .def_readwrite("prediction_error_control", &SmartEquilibriumTiming::prediction_error_control, "The time spent during on error control while searching during a smart prediction (in seconds).")
        .def_readwrite("prediction_taylor", &SmartEquilibriumTiming::prediction_taylor, "The time spent for the matrix-vector multiplication during a smart prediction (in seconds).")
        .def_readwrite("prediction_priority_update", &SmartEquilibriumTiming::prediction_priority_update, "The time spent for updating the priority related info of the clusters in the database (in seconds).")
        .def(py::self += py::self)
        ;

    py::class_<SmartEquilibriumResultDuringPrediction>(m, "SmartEquilibriumResultDuringPrediction")
        .def_readwrite("accepted", &SmartEquilibriumResultDuringPrediction::accepted)
        .def_readwrite("failed_with_species", &SmartEquilibriumResultDuringPrediction::failed_with_species)
        .def_readwrite("failed_with_amount", &SmartEquilibriumResultDuringPrediction::failed_with_amount)
        .def_readwrite("failed_with_chemical_potential", &SmartEquilibriumResultDuringPrediction::failed_with_chemical_potential)
        .def(py::self += py::self)
        ;

    py::class_<SmartEquilibriumResultDuringLearning>(m, "SmartEquilibriumResultDuringLearning")
        .def_readwrite("solve", &SmartEquilibriumResultDuringLearning::solve)
        .def(py::self += py::self)
        ;

    py::class_<SmartEquilibriumResult>(m, "SmartEquilibriumResult")
        .def("succeeded", &SmartEquilibriumResult::succeeded, "Return true if the calculation succeeded.")
        .def("predicted", &SmartEquilibriumResult::predicted, "Return true if the calculation was performed using a fast first-order Taylor prediction.")
        .def("learned", &SmartEquilibriumResult::learned, "Return true if the calculation was learned, not predicted, and performed using the conventional algorithm.")
        .def("failed", &SmartEquilibriumResult::failed, "Return true if the calculation failed.")
        .def("iterations", &SmartEquilibriumResult::iterations, "Return the number of iterations in the calculation.")
        .def_readwrite("prediction", &SmartEquilibriumResult::prediction)
        .def_readwrite("learning", &SmartEquilibriumResult::learning)
        .def_readwrite("timing", &SmartEquilibriumResult::timing)
        ;
}
