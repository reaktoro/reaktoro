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
#include <Reaktoro/Kinetics/KineticResult.hpp>

namespace Reaktoro {

void exportKineticResult(py::module& m)
{
    py::class_<KineticTiming>(m, "KineticTiming")
            .def_readwrite("solve", &KineticTiming::solve)
            .def_readwrite("initialize", &KineticTiming::initialize)
            .def_readwrite("integrate", &KineticTiming::integrate)
            .def_readwrite("integrate_chemical_properties", &KineticTiming::integrate_chemical_properties)
            .def_readwrite("integrate_reaction_rates", &KineticTiming::integrate_reaction_rates)
            .def_readwrite("integrate_sensitivity", &KineticTiming::integrate_sensitivity)
            .def_readwrite("integrate_equilibration", &KineticTiming::integrate_equilibration)
            .def_readwrite("equilibrate", &KineticTiming::equilibrate)
            .def(py::self += py::self)
            ;

    py::class_<KineticResult>(m, "KineticResult")
            .def_readwrite("timing", &KineticResult::timing)
            .def_readwrite("equilibrium", &KineticResult::equilibrium)
            .def_readwrite("smart_equilibrium", &KineticResult::smart_equilibrium)
            .def(py::self += py::self)
            ;
}

} // namespace Reaktoro
