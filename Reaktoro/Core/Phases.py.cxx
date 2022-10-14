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
        try { phases.add(phase.cast<GenericPhase const&>()); }
        catch(...)
        {
            try { phases.add(phase.cast<GenericPhasesGenerator const&>()); }
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

    m.def("speciate", speciate);

    py::class_<Exclude>(m, "Exclude")
        .def(py::init<>())
        .def_readwrite("tags", &Exclude::tags)
        ;

    m.def("exclude", exclude);

    py::class_<Phases>(m, "Phases")
        .def(py::init<const Database&>())
        .def(py::init(&rkt4py::createPhases))
        .def("add", py::overload_cast<const GenericPhase&>(&Phases::add))
        .def("add", py::overload_cast<const GenericPhasesGenerator&>(&Phases::add))
        .def("database", &Phases::database, return_internal_ref)
        .def("genericPhases", &Phases::genericPhases, return_internal_ref)
        .def("genericPhasesGenerators", &Phases::genericPhasesGenerators, return_internal_ref)
        .def("convert", &Phases::convert)
        ;

    py::implicitly_convertible<Phases, PhaseList>();

    py::class_<GenericPhase>(m, "GenericPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        .def("setName", &GenericPhase::setName, return_internal_ref)
        .def("setStateOfMatter", &GenericPhase::setStateOfMatter, return_internal_ref)
        .def("setAggregateState", &GenericPhase::setAggregateState, return_internal_ref)
        .def("setActivityModel", &GenericPhase::setActivityModel, return_internal_ref)
        .def("setIdealActivityModel", &GenericPhase::setIdealActivityModel, return_internal_ref)
        .def("named", &GenericPhase::named, return_internal_ref)
        .def("set", py::overload_cast<StateOfMatter>(&GenericPhase::set), return_internal_ref)
        .def("set", py::overload_cast<AggregateState>(&GenericPhase::set), return_internal_ref)
        .def("set", py::overload_cast<const ActivityModelGenerator&>(&GenericPhase::set), return_internal_ref)
        .def("name", &GenericPhase::name)
        .def("stateOfMatter", &GenericPhase::stateOfMatter)
        .def("aggregateState", &GenericPhase::aggregateState)
        .def("additionalAggregateStates", &GenericPhase::additionalAggregateStates, return_internal_ref)
        .def("species", &GenericPhase::species, return_internal_ref)
        .def("elements", &GenericPhase::elements, return_internal_ref)
        .def("activityModel", &GenericPhase::activityModel, return_internal_ref)
        .def("idealActivityModel", &GenericPhase::idealActivityModel, return_internal_ref)
        .def("convert", &GenericPhase::convert)
        ;

    py::class_<GenericPhasesGenerator>(m, "GenericPhasesGenerator")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        .def("setStateOfMatter", &GenericPhasesGenerator::setStateOfMatter, return_internal_ref)
        .def("setAggregateState", &GenericPhasesGenerator::setAggregateState, return_internal_ref)
        .def("setActivityModel", &GenericPhasesGenerator::setActivityModel, return_internal_ref)
        .def("setIdealActivityModel", &GenericPhasesGenerator::setIdealActivityModel, return_internal_ref)
        .def("set", py::overload_cast<StateOfMatter>(&GenericPhasesGenerator::set), return_internal_ref)
        .def("set", py::overload_cast<AggregateState>(&GenericPhasesGenerator::set), return_internal_ref)
        .def("set", py::overload_cast<const ActivityModelGenerator&>(&GenericPhasesGenerator::set), return_internal_ref)
        .def("stateOfMatter", &GenericPhasesGenerator::stateOfMatter)
        .def("aggregateState", &GenericPhasesGenerator::aggregateState)
        .def("additionalAggregateStates", &GenericPhasesGenerator::additionalAggregateStates, return_internal_ref)
        .def("species", &GenericPhasesGenerator::species, return_internal_ref)
        .def("elements", &GenericPhasesGenerator::elements, return_internal_ref)
        .def("activityModel", &GenericPhasesGenerator::activityModel, return_internal_ref)
        .def("idealActivityModel", &GenericPhasesGenerator::idealActivityModel, return_internal_ref)
        .def("convert", &GenericPhasesGenerator::convert)
        ;

    py::class_<AqueousPhase, GenericPhase>(m, "AqueousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<Speciate>())
        .def(py::init<Speciate, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<GaseousPhase, GenericPhase>(m, "GaseousPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<LiquidPhase, GenericPhase>(m, "LiquidPhase")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<SolidPhase, GenericPhase>(m, "SolidPhase")
        .def(py::init<const StringList&>())
        ;

    py::class_<MineralPhase, GenericPhase>(m, "MineralPhase")
        .def(py::init<String>())
        ;

    py::class_<MineralPhases, GenericPhasesGenerator>(m, "MineralPhases")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<CondensedPhase, GenericPhase>(m, "CondensedPhase")
        .def(py::init<String>())
        ;

    py::class_<CondensedPhases, GenericPhasesGenerator>(m, "CondensedPhases")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Speciate&>())
        .def(py::init<const Speciate&, const Exclude&>())
        .def(py::init<const Exclude&>())
        ;

    py::class_<IonExchangePhase, GenericPhase>(m, "IonExchangePhase")
        .def(py::init<String>())
        ;
}
