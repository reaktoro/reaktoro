//// Reaktoro is a unified framework for modeling chemically reactive systems.
////
//// Copyright (C) 2014-2015 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#pragma once
//
//// C++ includes
//#include <memory>
//
//// Reaktoro includes
//#include <Reaktoro/Common/ChemicalVector.hpp>
//#include <Reaktoro/Common/ThermoVector.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class Phase;
//
///// Defines a class for querying thermodynamic and chemical properties of a phase.
//class PhaseThermoProperties
//{
//public:
//    /// Construct a default PhaseThermoProperties instance.
//    PhaseThermoProperties();
//
//    /// Construct a PhaseThermoProperties instance with given Phase instance.
//    explicit PhaseThermoProperties(const Phase& phase);
//
//    /// Construct a copy of a PhaseThermoProperties instance.
//    PhaseThermoProperties(const PhaseThermoProperties& other);
//
//    /// Destroy this instance.
//    virtual ~PhaseThermoProperties();
//
//    /// Construct a copy of a PhaseThermoProperties instance.
//    auto operator=(PhaseThermoProperties other) -> PhaseThermoProperties&;
//
//    /// Update the thermodynamic properties of the phase.
//    /// @param T The new temperature (in units of K)
//    /// @param P The new pressure (in units of Pa)
//    auto update(double T, double P) -> void;
//
//    /// Return the temperature of the phase (in units of K).
//    auto temperature() const -> double;
//
//    /// Return the pressure of the phase (in units of Pa).
//    auto pressure() const -> double;
//
//    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
//    auto standardPartialMolarGibbsEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
//    auto standardPartialMolarEnthalpies() const -> ThermoVector;
//
//    /// Return the standard partial molar volumes of the species (in units of m3/mol).
//    auto standardPartialMolarVolumes() const -> ThermoVector;
//
//    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
//    auto standardPartialMolarEntropies() const -> ThermoVector;
//
//    /// Return the standard partial molar internal energies of the species (in units of J/mol).
//    auto standardPartialMolarInternalEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
//    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector;
//
//    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector;
//
//    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector;
//
//private:
//    struct Impl;
//
//    std::unique_ptr<Impl> pimpl;
//};
//
//} // namespace Reaktoro
