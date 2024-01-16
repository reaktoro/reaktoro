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
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Phases.hpp>
using namespace Reaktoro;

namespace rkt4py {

/// Creates a Phases object with given database and list of arguments expected to be Phase-like objects
auto createPhases(Database const& db, py::args gphases) -> Phases
{
    Phases phases(db);
    for(auto phase : gphases)
    {
        try { phases.add(phase.cast<GeneralPhase const&>()); }
        catch(...)
        {
            try { phases.add(phase.cast<GeneralPhasesGenerator const&>()); }
            catch(...)
            {
                errorif(true, "Could not create Phases with phase object:\n", py::str(phase));
            }
        }
    }
    return phases;
};

} // namespace rkt4py

void exportPhases(py::module& m)
{
    py::class_<Speciate>(m, "Speciate")
        .def(py::init<>())
        .def_readwrite("symbols", &Speciate::symbols)
        ;

    m.def("speciate", speciate, "The function used to specify phase species to be determined from element symbols in a list of substance formulas.");

    py::class_<Exclude>(m, "Exclude")
        .def(py::init<>())
        .def_readwrite("tags", &Exclude::tags)
        ;

    m.def("exclude", exclude, "The function used to specify species that should be filtered out when contructing a phase.");

    py::class_<Phases>(m, "Phases")
        .def(py::init<const Database&>())
        .def(py::init(&rkt4py::createPhases))
        .def("add", py::overload_cast<const GeneralPhase&>(&Phases::add))
        .def("add", py::overload_cast<const GeneralPhasesGenerator&>(&Phases::add))
        .def("database", &Phases::database, return_internal_ref)
        .def("generalPhases", &Phases::generalPhases, return_internal_ref)
        .def("generalPhasesGenerators", &Phases::generalPhasesGenerators, return_internal_ref)
        .def("convert", &Phases::convert)
        ;

    py::implicitly_convertible<Phases, PhaseList>();

    py::class_<GeneralPhase>(m, "GeneralPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        .def("setName", &GeneralPhase::setName, return_internal_ref)
        .def("setStateOfMatter", &GeneralPhase::setStateOfMatter, return_internal_ref)
        .def("setAggregateState", &GeneralPhase::setAggregateState, return_internal_ref)
        .def("setActivityModel", &GeneralPhase::setActivityModel, return_internal_ref)
        .def("setIdealActivityModel", &GeneralPhase::setIdealActivityModel, return_internal_ref)
        .def("named", &GeneralPhase::named, return_internal_ref)
        .def("set", py::overload_cast<StateOfMatter>(&GeneralPhase::set), return_internal_ref)
        .def("set", py::overload_cast<AggregateState>(&GeneralPhase::set), return_internal_ref)
        .def("set", py::overload_cast<const ActivityModelGenerator&>(&GeneralPhase::set), return_internal_ref)
        .def("name", &GeneralPhase::name)
        .def("stateOfMatter", &GeneralPhase::stateOfMatter)
        .def("aggregateState", &GeneralPhase::aggregateState)
        .def("additionalAggregateStates", &GeneralPhase::additionalAggregateStates, return_internal_ref)
        .def("species", &GeneralPhase::species, return_internal_ref)
        .def("elements", &GeneralPhase::elements, return_internal_ref)
        .def("activityModel", &GeneralPhase::activityModel, return_internal_ref)
        .def("idealActivityModel", &GeneralPhase::idealActivityModel, return_internal_ref)
        .def("convert", &GeneralPhase::convert)
        ;

    py::class_<GeneralPhasesGenerator>(m, "GeneralPhasesGenerator")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        .def("setStateOfMatter", &GeneralPhasesGenerator::setStateOfMatter, return_internal_ref)
        .def("setAggregateState", &GeneralPhasesGenerator::setAggregateState, return_internal_ref)
        .def("setActivityModel", &GeneralPhasesGenerator::setActivityModel, return_internal_ref)
        .def("setIdealActivityModel", &GeneralPhasesGenerator::setIdealActivityModel, return_internal_ref)
        .def("set", py::overload_cast<StateOfMatter>(&GeneralPhasesGenerator::set), return_internal_ref)
        .def("set", py::overload_cast<AggregateState>(&GeneralPhasesGenerator::set), return_internal_ref)
        .def("set", py::overload_cast<const ActivityModelGenerator&>(&GeneralPhasesGenerator::set), return_internal_ref)
        .def("stateOfMatter", &GeneralPhasesGenerator::stateOfMatter)
        .def("aggregateState", &GeneralPhasesGenerator::aggregateState)
        .def("additionalAggregateStates", &GeneralPhasesGenerator::additionalAggregateStates, return_internal_ref)
        .def("species", &GeneralPhasesGenerator::species, return_internal_ref)
        .def("elements", &GeneralPhasesGenerator::elements, return_internal_ref)
        .def("activityModel", &GeneralPhasesGenerator::activityModel, return_internal_ref)
        .def("idealActivityModel", &GeneralPhasesGenerator::idealActivityModel, return_internal_ref)
        .def("convert", &GeneralPhasesGenerator::convert)
        ;

    py::class_<AqueousPhase, GeneralPhase>(m, "AqueousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<Speciate>())
        .def(py::init<Speciate, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<GaseousPhase, GeneralPhase>(m, "GaseousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<LiquidPhase, GeneralPhase>(m, "LiquidPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<SolidPhase, GeneralPhase>(m, "SolidPhase")
        .def(py::init<const StringList&>())
        ;

    py::class_<MineralPhase, GeneralPhase>(m, "MineralPhase")
        .def(py::init<String>())
        ;

    py::class_<MineralPhases, GeneralPhasesGenerator>(m, "MineralPhases")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<CondensedPhase, GeneralPhase>(m, "CondensedPhase")
        .def(py::init<String>())
        ;

    py::class_<CondensedPhases, GeneralPhasesGenerator>(m, "CondensedPhases")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<IonExchangePhase, GeneralPhase>(m, "IonExchangePhase")
        .def(py::init<String>())
        ;
}
