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
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Math/ODE.hpp>

namespace Reaktoro {

void exportODE(py::module& m)
{
    py::enum_<ODEStepMode>(m, "ODEStepMode")
        .value("Adams", ODEStepMode::Adams)
        .value("BDF", ODEStepMode::BDF)
        ;

    py::enum_<ODEIterationMode>(m, "ODEIterationMode")
        .value("Functional", ODEIterationMode::Functional)
        .value("Newton", ODEIterationMode::Newton)
        ;

    py::class_<ODEOptions>(m, "ODEOptions")
        .def(py::init<>())
        .def_readwrite("step", &ODEOptions::step)
        .def_readwrite("stability_limit_detection", &ODEOptions::stability_limit_detection)
        .def_readwrite("initial_step", &ODEOptions::initial_step)
        .def_readwrite("stop_time", &ODEOptions::stop_time)
        .def_readwrite("min_step", &ODEOptions::min_step)
        .def_readwrite("max_step", &ODEOptions::max_step)
        .def_readwrite("reltol", &ODEOptions::reltol)
        .def_readwrite("abstol", &ODEOptions::abstol)
        .def_readwrite("max_order_bdf", &ODEOptions::max_order_bdf)
        .def_readwrite("max_order_adams", &ODEOptions::max_order_adams)
        .def_readwrite("max_num_steps", &ODEOptions::max_num_steps)
        .def_readwrite("max_hnil_warnings", &ODEOptions::max_hnil_warnings)
        .def_readwrite("max_num_error_test_failures", &ODEOptions::max_num_error_test_failures)
        .def_readwrite("max_num_nonlinear_iterations", &ODEOptions::max_num_nonlinear_iterations)
        .def_readwrite("max_num_convergence_failures", &ODEOptions::max_num_convergence_failures)
        .def_readwrite("nonlinear_convergence_coefficient", &ODEOptions::nonlinear_convergence_coefficient)
        .def_readwrite("abstols", &ODEOptions::abstols)
        ;

//
//    auto function1 = static_cast<const ODEFunction&(ODEProblem::*)() const>(&ODEProblem::function);
//    auto function2 = static_cast<int(ODEProblem::*)(double, VectorConstRef, VectorReff) const>(&ODEProblem::function);
//
//    auto jacobian1 = static_cast<const ODEJacobian&(ODEProblem::*)() const>(&ODEProblem::jacobian);
//    auto jacobian2 = static_cast<int(ODEProblem::*)(double, VectorConstRef, MatrixReff) const>(&ODEProblem::jacobian);
//
//    py::class_<ODEProblem>(m, "ODEProblem")
//        .def(py::init<>())
//        .def("setNumEquations", &ODEProblem::setNumEquations)
//        .def("setFunction", &ODEProblem::setFunction)
//        .def("setJacobian", &ODEProblem::setJacobian)
//        .def("initialized", &ODEProblem::initialized)
//        .def("numEquations", &ODEProblem::numEquations)
//        .def("function", function1, py::return_value_policy::reference_internal)
//        .def("jacobian", jacobian1, py::return_value_policy::reference_internal)
//        .def("function", function2)
//        .def("jacobian", jacobian2)
//        ;
//
//    auto integrate1 = static_cast<void(ODESolver::*)(double&, VectorRef)>(&ODESolver::integrate);
//    auto integrate2 = static_cast<void(ODESolver::*)(double&, VectorRef, double)>(&ODESolver::integrate);
//
//    py::class_<ODESolver>(m, "ODESolver")
//        .def(py::init<>())
//        .def("setOptions", &ODESolver::setOptions)
//        .def("setProblem", &ODESolver::setProblem)
//        .def("initialize", &ODESolver::initialize)
//        .def("integrate", integrate1)
//        .def("integrate", integrate2)
//        .def("solve", &ODESolver::solve)
//        ;
}

} // namespace Reaktoro
