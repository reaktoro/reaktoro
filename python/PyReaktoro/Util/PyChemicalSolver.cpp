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

#include "PyChemicalSolver.hpp"

// Numpy includes
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>

// Boost includes
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Util/ChemicalField.hpp>
#include <Reaktoro/Util/ChemicalSolver.hpp>

namespace Reaktoro {
namespace PyChemicalSolver {

auto setStatesHelper1(ChemicalSolver& self, const py::object& states) -> bool
{
    py::extract<const ChemicalState&> state(states);
    if(state.check())
    {
        self.setStates(state());
        return true;
    }
    return false;
}

auto setStatesHelper2(ChemicalSolver& self, const py::object& states) -> bool
{
    const Index len = py::len(states);
    Assert(len == self.numPoints(),
        "Could not set the chemical states at every field point.",
        "Expecting the same number of chemical states as there are field points.");
    for(Index i = 0; i < len; ++i)
    {
        py::extract<const ChemicalState&> state(states[i]);
        Assert(state.check(),
            "Could not set the chemical states with given chemical states.",
            "Expecting chemical states with type reaktoro.ChemicalState.");
        self.setStateAt(i, state());
    }
    return true;
}

auto setStates(ChemicalSolver& self, const py::object& states) -> void
{
    if(setStatesHelper1(self, states)) return;
    if(setStatesHelper2(self, states)) return;
}

auto setStateAtHelper1(ChemicalSolver& self, const py::object& ipoints, const py::object& states) -> bool
{
    py::extract<int> index(ipoints);
    py::extract<const ChemicalState&> state(states);
    if(index.check() && state.check())
    {
        self.setStateAt(index(), state());
        return true;
    }
    return false;
}

auto setStateAtHelper2(ChemicalSolver& self, const py::object& ipoints, const py::object& states) -> bool
{
    const Index len = py::len(ipoints);
    Assert(len <= self.numPoints(),
        "Could not set the chemical state at given field points.",
        "Expecting number of indices not greater than the number of field points.");
    py::extract<const ChemicalState&> state(states);
    if(state.check())
    {
         for(Index i = 0; i < len; ++i)
         {
             py::extract<int> index(ipoints[i]);
             Assert(index.check(),
                 "Could not set the chemical states with the given field point indices.",
                 "Expecting indices of integer type.");
             self.setStateAt(index(), state());
         }
        return true;
    }
    return false;
}

auto setStateAtHelper3(ChemicalSolver& self, const py::object& ipoints, const py::object& states) -> bool
{
    const Index len_ipoints = py::len(ipoints);
    const Index len_states  = py::len(states);

    Assert(len_states == len_ipoints,
        "Could not set the chemical state at given field points.",
        "Expecting the same number of field point indices and chemical states.");

    Assert(len_ipoints <= self.numPoints(),
        "Could not set the chemical state at given field points.",
        "Expecting number of indices not greater that the number of field points.");

    for(Index i = 0; i < len_ipoints; ++i)
    {
        py::extract<int> index(ipoints[i]);
        py::extract<const ChemicalState&> state(states[i]);
        Assert(index.check(),
            "Could not set the chemical states with the given field point indices.",
            "Expecting indices of integer type.");
        Assert(state.check(),
            "Could not set the chemical states with given chemical states.",
            "Expecting chemical states with type reaktoro.ChemicalState.");
        self.setStateAt(index(), state());
    }
    return true;
}

auto setStateAt(ChemicalSolver& self, const py::object& ipoints, const py::object& states) -> void
{
    if(setStateAtHelper1(self, ipoints, states)) return;
    if(setStateAtHelper2(self, ipoints, states)) return;
    if(setStateAtHelper3(self, ipoints, states)) return;
}

auto equilibrate(
        ChemicalSolver& self,
        const py::numeric::array& T,
        const py::numeric::array& P,
        const py::numeric::array& b) -> void
{
    PyArrayObject* ptr_T = (PyArrayObject*)T.ptr();
    PyArrayObject* ptr_P = (PyArrayObject*)P.ptr();
    PyArrayObject* ptr_b = (PyArrayObject*)b.ptr();

    Assert(PyArray_ISFLOAT(ptr_T),
        "Could not perform equilibrium calculations.",
        "Expecting a numpy.array with dtype float64 for temperatures.");

    Assert(PyArray_ISFLOAT(ptr_P),
        "Could not perform equilibrium calculations.",
        "Expecting a numpy.array with dtype float64 for pressures.");

    Assert(PyArray_ISFLOAT(ptr_b),
        "Could not perform equilibrium calculations.",
        "Expecting a numpy.array with dtype float64 for element amounts.");

    const Index len_T = py::len(T);
    const Index len_P = py::len(P);
    const Index len_b = py::len(b);

    double* data_T = (double*)PyArray_DATA(ptr_T);
    double* data_P = (double*)PyArray_DATA(ptr_P);
    double* data_b = (double*)PyArray_DATA(ptr_b);

    ChemicalSolver::Array<double> array_T(data_T, len_T);
    ChemicalSolver::Array<double> array_P(data_P, len_P);
    ChemicalSolver::Array<double> array_b(data_b, len_b);

    self.equilibrate(array_T, array_P, array_b);
}

} // namespace

auto export_ChemicalSolver() -> void
{
    import_array();

    py::class_<ChemicalSolver>("ChemicalSolver")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&, Index>())
        .def(py::init<const ReactionSystem&, Index>())
        .def("numPoints", &ChemicalSolver::numPoints)
        .def("numEquilibriumElements", &ChemicalSolver::numEquilibriumElements)
        .def("numKineticSpecies", &ChemicalSolver::numKineticSpecies)
        .def("numComponents", &ChemicalSolver::numComponents)
        .def("setPartition", &ChemicalSolver::setPartition)
        .def("setStates", PyChemicalSolver::setStates)
        .def("setStateAt", PyChemicalSolver::setStateAt)
        .def("equilibrate", PyChemicalSolver::equilibrate)
        .def("react", &ChemicalSolver::react)
        .def("state", &ChemicalSolver::state, py::return_internal_reference<>())
        .def("states", &ChemicalSolver::states, py::return_internal_reference<>())
        .def("componentAmounts", &ChemicalSolver::componentAmounts, py::return_internal_reference<>())
        .def("equilibriumSpeciesAmounts", &ChemicalSolver::equilibriumSpeciesAmounts, py::return_internal_reference<>())
        .def("porosity", &ChemicalSolver::porosity, py::return_internal_reference<>())
        .def("fluidSaturations", &ChemicalSolver::fluidSaturations, py::return_internal_reference<>())
        .def("fluidDensities", &ChemicalSolver::fluidDensities, py::return_internal_reference<>())
        .def("fluidVolumes", &ChemicalSolver::fluidVolumes, py::return_internal_reference<>())
        .def("fluidTotalVolume", &ChemicalSolver::fluidTotalVolume, py::return_internal_reference<>())
        .def("solidTotalVolume", &ChemicalSolver::solidTotalVolume, py::return_internal_reference<>())
        .def("componentRates", &ChemicalSolver::componentRates, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
