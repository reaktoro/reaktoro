// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
using namespace Reaktoro;

void exportChemicalPropsPhase(py::module& m)
{
    auto update1 = [](ChemicalPropsPhase& self, const real& T, const real& P, ArrayXrConstRef n, Map<String, Any>& extra)
    {
        self.update(T, P, n, extra);
    };

    auto updateIdeal1 = [](ChemicalPropsPhase& self, const real& T, const real& P, ArrayXrConstRef n, Map<String, Any>& extra)
    {
        self.updateIdeal(T, P, n, extra);
    };

    py::class_<ChemicalPropsPhase>(m, "ChemicalPropsPhase")
        .def(py::init<const Phase&>())
        .def("update", update1, "Update the chemical properties of the phase.")
        .def("updateIdeal", updateIdeal1, "Update the chemical properties of the phase using ideal activity models.")
        .def("phase", &ChemicalPropsPhase::phase, py::return_value_policy::reference_internal, "Return the underlying Phase object.")
        .def("data", &ChemicalPropsPhase::data, py::return_value_policy::reference_internal, "Return the primary chemical property data of the phase from which others are calculated.")
        .def("temperature", &ChemicalPropsPhase::temperature, "Return the temperature of the phase (in K).")
        .def("pressure", &ChemicalPropsPhase::pressure, "Return the pressure of the phase (in Pa).")
        .def("speciesAmounts", &ChemicalPropsPhase::speciesAmounts, py::return_value_policy::reference_internal, "Return the amounts of the species in the phase (in mol).")
        .def("moleFractions", &ChemicalPropsPhase::moleFractions, py::return_value_policy::reference_internal, "Return the mole fractions of the species in the phase.")
        .def("lnActivityCoefficients", &ChemicalPropsPhase::lnActivityCoefficients, py::return_value_policy::reference_internal, "Return the ln activity coefficients of the species in the phase.")
        .def("lnActivities", &ChemicalPropsPhase::lnActivities, py::return_value_policy::reference_internal, "Return the ln activities of the species in the phase.")
        .def("chemicalPotentials", &ChemicalPropsPhase::chemicalPotentials, py::return_value_policy::reference_internal, "Return the chemical potentials of the species in the phase (in J/mol).")
        .def("standardVolumes", &ChemicalPropsPhase::standardVolumes, py::return_value_policy::reference_internal, "Return the standard partial molar volumes of the species in the phase (in m³/mol).")
        .def("standardVolumesT", &ChemicalPropsPhase::standardVolumesT, py::return_value_policy::reference_internal, "Return the temperature derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).")
        .def("standardVolumesP", &ChemicalPropsPhase::standardVolumesP, py::return_value_policy::reference_internal, "Return the pressure derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).")
        .def("standardGibbsEnergies", &ChemicalPropsPhase::standardGibbsEnergies, py::return_value_policy::reference_internal, "Return the standard partial molar Gibbs energies of formation of the species in the phase (in J/mol).")
        .def("standardEnthalpies", &ChemicalPropsPhase::standardEnthalpies, py::return_value_policy::reference_internal, "Return the standard partial molar enthalpies of formation of the species in the phase (in J/mol).")
        .def("standardEntropies", &ChemicalPropsPhase::standardEntropies, "Return the standard partial molar entropies of formation of the species in the phase (in J/(mol·K)).")
        .def("standardInternalEnergies", &ChemicalPropsPhase::standardInternalEnergies, "Return the standard partial molar internal energies of formation of the species in the phase (in J/mol).")
        .def("standardHelmholtzEnergies", &ChemicalPropsPhase::standardHelmholtzEnergies, "Return the standard partial molar Helmholtz energies of formation of the species in the phase (in J/mol).")
        .def("standardHeatCapacitiesConstP", &ChemicalPropsPhase::standardHeatCapacitiesConstP, py::return_value_policy::reference_internal, "Return the standard partial molar isobaric heat capacities of the species in the phase (in J/(mol·K)).")
        .def("standardHeatCapacitiesConstV", &ChemicalPropsPhase::standardHeatCapacitiesConstV, py::return_value_policy::reference_internal, "Return the standard partial molar isochoric heat capacities of the species in the phase (in J/(mol·K)).")
        .def("molarMass", &ChemicalPropsPhase::molarMass, "Return the molar mass of the phase (in kg/mol).")
        .def("molarVolume", &ChemicalPropsPhase::molarVolume, "Return the molar volume of the phase (in m³/mol).")
        .def("molarVolumeT", &ChemicalPropsPhase::molarVolumeT, "Return the temperature derivative of the molar volume of the phase (in m³/(mol·K)).")
        .def("molarVolumeP", &ChemicalPropsPhase::molarVolumeP, "Return the pressure derivative of the molar volume of the phase (in m³/(mol·Pa)).")
        .def("molarGibbsEnergy", &ChemicalPropsPhase::molarGibbsEnergy, "Return the molar Gibbs energy of formation of the phase (in J/mol).")
        .def("molarEnthalpy", &ChemicalPropsPhase::molarEnthalpy, "Return the molar enthalpy of formation of the phase (in J/mol).")
        .def("molarEntropy", &ChemicalPropsPhase::molarEntropy, "Return the molar entropy of formation of the phase (in J/(mol·K)).")
        .def("molarInternalEnergy", &ChemicalPropsPhase::molarInternalEnergy, "Return the molar internal energy of formation of the phase (in J/mol).")
        .def("molarHelmholtzEnergy", &ChemicalPropsPhase::molarHelmholtzEnergy, "Return the molar Helmholtz energy of formation of the phase (in J/mol).")
        .def("molarHeatCapacityConstP", &ChemicalPropsPhase::molarHeatCapacityConstP, "Return the molar isobaric heat capacity of the phase (in J/(mol·K)).")
        .def("molarHeatCapacityConstV", &ChemicalPropsPhase::molarHeatCapacityConstV, "Return the molar isochoric heat capacity of the phase (in J/(mol·K)).")
        .def("density", &ChemicalPropsPhase::density, "Return the density of the phase (in kg/m³).")
        .def("amount", &ChemicalPropsPhase::amount, "Return the sum of species amounts in the phase (in mol).")
        .def("mass", &ChemicalPropsPhase::mass, "Return the sum of species masses in the phase (in kg).")
        .def("volume", &ChemicalPropsPhase::volume, "Return the volume of the phase (in m³).")
        .def("volumeT", &ChemicalPropsPhase::volumeT, "Return the temperature derivative of the volume of the phase (in m³/K).")
        .def("volumeP", &ChemicalPropsPhase::volumeP, "Return the pressure derivative of the volume of the phase (in m³/Pa).")
        .def("gibbsEnergy", &ChemicalPropsPhase::gibbsEnergy, "Return the Gibbs energy of the phase (in J).")
        .def("enthalpy", &ChemicalPropsPhase::enthalpy, "Return the enthalpy of the phase (in J).")
        .def("entropy", &ChemicalPropsPhase::entropy, "Return the entropy of the phase (in J/K).")
        .def("internalEnergy", &ChemicalPropsPhase::internalEnergy, "Return the internal energy of the phase (in J).")
        .def("helmholtzEnergy", &ChemicalPropsPhase::helmholtzEnergy, "Return the Helmholtz energy of the phase (in J).")
        .def("heatCapacityConstP", &ChemicalPropsPhase::heatCapacityConstP, "Return the isobaric heat capacity of the phase (in J/K).")
        .def("heatCapacityConstV", &ChemicalPropsPhase::heatCapacityConstV, "Return the isochoric heat capacity of the phase (in J/K).")
        .def("soundSpeed", &ChemicalPropsPhase::soundSpeed, "Return the speed of sound in the phase (in m/s).")
        ;
}
