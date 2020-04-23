// Reaktoro is a unified framework for modeling chemically reactive systems.
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
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseSolver.hpp>

namespace Reaktoro {

    void exportEquilibriumInverseSolver(py::module& m)
    {
        py::class_<EquilibriumInverseSolver>(m, "EquilibriumInverseSolver")
                .def(py::init<const ChemicalSystem&>())
                .def(py::init<const Partition&>())
                .def("setOptions", &EquilibriumInverseSolver::setOptions)
                .def("solve", &EquilibriumInverseSolver::solve)
                .def("sensitivity", &EquilibriumInverseSolver::sensitivity, py::return_value_policy::reference_internal)
                ;
    }

} // namespace Reaktoro
