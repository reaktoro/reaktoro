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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ScalarTypes.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
struct SpeciesElectroState;
struct SpeciesThermoState;
struct WaterElectroState;
struct WaterThermoState;

/// Calculate the thermodynamic state of solvent water using the HKF model.
auto speciesThermoStateSolventHKF(Temperature T, Pressure P, const WaterThermoState& wts) -> SpeciesThermoState;

/// Calculate the thermodynamic state of an aqueous solute using the HKF model.
auto speciesThermoStateSoluteHKF(Temperature T, Pressure P, const AqueousSpecies& species, const SpeciesElectroState& aes, const WaterElectroState& wes) -> SpeciesThermoState;

/// Calculate the thermodynamic state of an aqueous species using the HKF model.
auto speciesThermoStateHKF(Temperature T, Pressure P, const AqueousSpecies& species) -> SpeciesThermoState;

/// Calculate the thermodynamic state of a gaseous species using the HKF model.
auto speciesThermoStateHKF(Temperature T, Pressure P, const GaseousSpecies& species) -> SpeciesThermoState;

/// Calculate the thermodynamic state of a mineral species using the HKF model.
auto speciesThermoStateHKF(Temperature T, Pressure P, const MineralSpecies& species) -> SpeciesThermoState;

} // namespace Reaktoro
