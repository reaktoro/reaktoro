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

#pragma once

namespace Reaktoro {

// Forward declarations
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
class Pressure;
class Temperature;
class ThermoScalar;
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
