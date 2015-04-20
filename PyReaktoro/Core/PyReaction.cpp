// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PyReaction.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>

// PyReator includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

auto export_Reaction() -> void
{
    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    py::class_<Reaction>("Reaction")
        .def(py::init<>())
        .def("setName", &Reaction::setName)
        .def("setEquilibriumConstantFunction", &Reaction::setEquilibriumConstantFunction)
        .def("setStandardGibbsEnergyFunction", &Reaction::setStandardGibbsEnergyFunction)
        .def("setStandardHelmholtzEnergyFunction", &Reaction::setStandardHelmholtzEnergyFunction)
        .def("setStandardInternalEnergyFunction", &Reaction::setStandardInternalEnergyFunction)
        .def("setStandardEnthalpyFunction", &Reaction::setStandardEnthalpyFunction)
        .def("setStandardEntropyFunction", &Reaction::setStandardEntropyFunction)
        .def("setStandardVolumeFunction", &Reaction::setStandardVolumeFunction)
        .def("setStandardHeatCapacityFunction", &Reaction::setStandardHeatCapacityFunction)
        .def("setRate", &Reaction::setRate)
        .def("name", &Reaction::name)
        .def("equilibriumConstantFunction", &Reaction::equilibriumConstantFunction, return_const_ref())
        .def("standardGibbsEnergyFunction", &Reaction::standardGibbsEnergyFunction, return_const_ref())
        .def("standardHelmholtzEnergyFunction", &Reaction::standardHelmholtzEnergyFunction, return_const_ref())
        .def("standardInternalEnergyFunction", &Reaction::standardInternalEnergyFunction, return_const_ref())
        .def("standardEnthalpyFunction", &Reaction::standardEnthalpyFunction, return_const_ref())
        .def("standardEntropyFunction", &Reaction::standardEntropyFunction, return_const_ref())
        .def("standardVolumeFunction", &Reaction::standardVolumeFunction, return_const_ref())
        .def("standardHeatCapacityFunction", &Reaction::standardHeatCapacityFunction, return_const_ref())
        .def("rateFunction", &Reaction::rateFunction, return_const_ref())
        .def("equation", &Reaction::equation, return_const_ref())
        .def("system", &Reaction::system, return_const_ref())
        .def("species", &Reaction::species, return_const_ref())
        .def("indices", &Reaction::indices, return_const_ref())
        .def("stoichiometries", &Reaction::stoichiometries, return_const_ref())
        .def("stoichiometry", &Reaction::stoichiometry)
        .def("lnEquilibriumConstant", &Reaction::lnEquilibriumConstant)
        .def("standardGibbsEnergy", &Reaction::standardGibbsEnergy)
        .def("standardHelmholtzEnergy", &Reaction::standardHelmholtzEnergy)
        .def("standardInternalEnergy", &Reaction::standardInternalEnergy)
        .def("standardEnthalpy", &Reaction::standardEnthalpy)
        .def("standardEntropy", &Reaction::standardEntropy)
        .def("standardVolume", &Reaction::standardVolume)
        .def("standardHeatCapacity", &Reaction::standardHeatCapacity)
        .def("rate", &Reaction::rate)
        .def("lnReactionQuotient", &Reaction::lnReactionQuotient)
        ;

    export_std_vector<Reaction>("ReactionVector");
}

} // namespace Reaktoro
