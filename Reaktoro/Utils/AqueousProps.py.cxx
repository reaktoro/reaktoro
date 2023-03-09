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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>
using namespace Reaktoro;

void exportAqueousProps(py::module& m)
{
    auto saturationIndexLn   = [](AqueousProps const& self, StringOrIndex const& species) -> real { errorif(true, "Method AqueousProps::saturationIndexLn has been deprecated. Rely on the use of saturationIndex(species) instead."); return {}; };
    auto saturationIndexLg   = [](AqueousProps const& self, StringOrIndex const& species) -> real { errorif(true, "Method AqueousProps::saturationIndexLg has been deprecated. Rely on the use of saturationIndex(species) instead."); return {}; };
    auto saturationIndicesLn = [](AqueousProps const& self) -> real { errorif(true, "Method AqueousProps::saturationIndicesLn has been deprecated. Rely on the use of saturationIndices() instead."); return {}; };
    auto saturationIndicesLg = [](AqueousProps const& self) -> real { errorif(true, "Method AqueousProps::saturationIndicesLg has been deprecated. Rely on the use of saturationIndices() instead."); return {}; };

    py::class_<AqueousProps>(m, "AqueousProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def(py::init<const ChemicalProps&>())
        .def_static("compute", &AqueousProps::compute, "Compute an AqueousProps object with given ChemicalProps object.")
        .def("setActivityModel", &AqueousProps::setActivityModel, "Set an activity model for a non-aqueous species that will be used in the calculation of its saturation index.")
        .def("update", py::overload_cast<const ChemicalState&>(&AqueousProps::update), "Update the aqueous properties with given chemical state of the system.")
        .def("update", py::overload_cast<const ChemicalProps&>(&AqueousProps::update), "Update the aqueous properties with given chemical properties of the system.")
        .def("temperature", &AqueousProps::temperature, "Return the temperature of the aqueous phase (in K).")
        .def("pressure", &AqueousProps::pressure, "Return the pressure of the aqueous phase (in Pa).")
        .def("waterAmount", &AqueousProps::waterAmount, "Return the amount of solvent water in the aqueous phase (in mol).")
        .def("waterMass", &AqueousProps::waterMass, "Return the mass of solvent water in the aqueous phase (in kg).")
        .def("charge", &AqueousProps::charge, "Return the electric charge in the aqueous phase (in mol).")
        .def("chargeMolality", &AqueousProps::chargeMolality, "Return the molality concentration in the aqueous phase (in molal).")
        .def("elementMolality", &AqueousProps::elementMolality, "Return the molality of an element (in molal).")
        .def("elementMolalities", &AqueousProps::elementMolalities, "Return the molality concentrations of the elements in  (in molal).")
        .def("speciesMolality", &AqueousProps::speciesMolality, "Return the molality of an aqueous solute species (in molal).")
        .def("speciesMolalities", &AqueousProps::speciesMolalities, "Return the molality concentrations of the species (in molal).")
        .def("ionicStrength", &AqueousProps::ionicStrength, "Return the effective ionic strength of the aqueous phase (in molal). Equivalent to @ref ionicStrengthEffective.")
        .def("ionicStrengthEffective", &AqueousProps::ionicStrengthEffective, "Return the effective ionic strength of the aqueous phase (in molal).")
        .def("ionicStrengthStoichiometric", &AqueousProps::ionicStrengthStoichiometric, "Return the stoichiometric ionic strength of the aqueous phase (in molal).")
        .def("pH", &AqueousProps::pH, "Return the pH of the aqueous phase.")
        .def("pE", &AqueousProps::pE, "Return the pE of the aqueous phase.")
        .def("Eh", &AqueousProps::Eh, "Return the reduction potential of the aqueous phase (in V).")
        .def("alkalinity", &AqueousProps::alkalinity, "Return the total alkalinity of the aqueous phase (in eq/L).")
        .def("saturationSpecies", &AqueousProps::saturationSpecies, "Return the non-aqueous species that could be formed from the aqueous solution.")
        .def("saturationIndex", &AqueousProps::saturationIndex, "Return the saturation index SI ≡ log(Ω) = log(IAP/K) of a non-aqueous species.")
        .def("saturationIndices", &AqueousProps::saturationIndices, "Return the saturation indices of all non-aqueous species.")
        .def("saturationRatio", &AqueousProps::saturationRatio, "Return the saturation ratio SR ≡ Ω = IAP/K of a non-aqueous species.")
        .def("saturationRatios", &AqueousProps::saturationRatios, "Return the saturation ratios of all non-aqueous species.")
        .def("saturationRatiosLn", &AqueousProps::saturationRatiosLn, "Return the saturation ratios of all non-aqueous species (in natural log).")
        .def("props", &AqueousProps::props, return_internal_ref, "Return the underlying ChemicalProps object.")
        .def("system", &AqueousProps::system, return_internal_ref, "Return the underlying ChemicalSystem object.")
        .def("phase", &AqueousProps::phase, return_internal_ref, "Return the underlying Phase object for the aqueous phase.")
        .def("output", py::overload_cast<std::ostream&>(&AqueousProps::output, py::const_), "Output the properties of the aqueous phase to a stream.")
        .def("output", py::overload_cast<const String&>(&AqueousProps::output, py::const_), "Output the properties of the aqueous phase to a file.")
        .def("__repr__", [](const AqueousProps& self) { std::stringstream ss; ss << self; return ss.str(); })

        // DEPRECATED METHODS TO BE REMOVED IN THE NEAR FUTURE

        .def("saturationIndexLn", saturationIndexLn, "Return the saturation index of a given species (in natural log).")
        .def("saturationIndexLg", saturationIndexLg, "Return the saturation index of a given species (in log base 10).")
        .def("saturationIndicesLn", saturationIndicesLn, "Return the saturation indices of all non-aqueous species (in natural log).")
        .def("saturationIndicesLg", saturationIndicesLg, "Return the saturation indices of all non-aqueous species (in log base 10).")
        ;
}
