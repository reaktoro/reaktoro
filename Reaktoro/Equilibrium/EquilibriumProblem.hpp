// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// #pragma once

// // Reaktoro includes
// #include <Reaktoro/Common/Types.hpp>
// #include <Reaktoro/Common/Matrix.hpp>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalState;
// class ChemicalSystem;

// /// The class used to define an equilibrium problem.
// class EquilibriumProblem
// {
// public:

//     // Forward declarations
//     class Control;
//     class Until;


//     /// Construct an EquilibriumProblem instance
//     explicit EquilibriumProblem(const ChemicalSystem& system);

//     /// Construct a copy of a EquilibriumProblem instance
//     EquilibriumProblem(const EquilibriumProblem& other);

//     /// Destroy this EquilibriumProblem instance
//     virtual ~EquilibriumProblem();

//     /// Assign an EquilibriumProblem instance to this
//     auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

//     /// Set the temperature for the equilibrium calculation.
//     /// By default, the temperature is 25 &deg;C.
//     /// @param val The temperature value (in unit of K)
//     auto setTemperature(double val) -> EquilibriumProblem&;

//     /// Set the temperature for the equilibrium calculation with given unit.
//     /// By default, the temperature is 25 &deg;C.
//     /// @param val The temperature value
//     /// @param unit The unit of the temperature (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine)
//     auto setTemperature(double val, std::string unit) -> EquilibriumProblem&;

//     /// Set the pressure for the equilibrium calculation.
//     /// By default, the pressure is 1 bar.
//     /// @param val The pressure value (in unit of Pa).
//     auto setPressure(double val) -> EquilibriumProblem&;

//     /// Set the pressure for the equilibrium calculation.
//     /// By default, the pressure is 1 bar.
//     /// @param val The pressure value
//     /// @param unit The unit of the pressure (Pa, kPa, MPa, GPa, atm, mmHg, inHg, psi, kpsi, Mpsi, psf, bar, torr, inH2O, ftH2O, pascal)
//     auto setPressure(double val, std::string unit) -> EquilibriumProblem&;

//     /// Set the mole amounts of each element for the equilibrium calculation.
//     /// @param b The vector of mole amounts of each element (in unit of mol)
//     auto setElementAmounts(VectorXrConstRef b) -> EquilibriumProblem&;

//     /// Set the mole amounts of each element for the equilibrium calculation.
//     /// @param amount The mole amount for all elements (in unit of mol)
//     auto setElementAmounts(double amount) -> EquilibriumProblem&;

//     /// Set the mole amount of an element for the equilibrium calculation (in unit of mol)
//     /// @param ielement The index of the element
//     /// @param amount The same mole amount for all elements (in unit of mol)
//     auto setElementAmount(Index ielement, double amount) -> EquilibriumProblem&;

//     /// Set the mole amount of an element for the equilibrium calculation (in unit of mol)
//     /// @param element The name of the element
//     /// @param amount The same mole amount for all elements (in unit of mol)
//     auto setElementAmount(std::string element, double amount) -> EquilibriumProblem&;

//     /// Set the mole amount of electrical charge.
//     /// @param amount The mole amount of electrical charge (in unit of mol)
//     auto setElectricalCharge(double amount) -> EquilibriumProblem&;

//     /// Add a given amount of a compound or species to the equilibrium recipe.
//     /// This method will first check if the given compound is present in the chemical system.
//     /// If true, then `addSpecies` will be called. Otherwise, `addCompound` will be called.
//     /// @param name The name of the compound or species
//     /// @param amount The amount of the compound or species
//     /// @param unit The unit of the amount (must be convertible to either mol or kg)
//     auto add(std::string name, double amount, std::string unit) -> EquilibriumProblem&;

//     /// Add the mole amounts of the equilibrium species in a ChemicalState instance to the equilibrium recipe.
//     /// This method only extracts the mole amounts of equilibrium species in the given chemical state.
//     /// If not all species in the system is under equilibrium assumption for this equilibrium calculation,
//     /// ensure you set the partition of the chemical system before calling this method.
//     /// @note If a multiplication factor is needed, for example `2.0`, use `add(2.0*state)`.
//     /// @param state The ChemicalState instance with the mole amounts of the species
//     /// @see setPartition
//     auto add(const ChemicalState& state) -> EquilibriumProblem&;

//     /// Add a given amount of a compound to the equilibrium recipe.
//     /// The compound must not have a chemical element that is not present in the chemical system.
//     /// @param name The name of the compound (e.g., H2O, CaCO3)
//     /// @param amount The amount of the compound
//     /// @param unit The unit of the amount (must be convertible to either mol or kg)
//     auto addCompound(std::string name, double amount, std::string unit) -> EquilibriumProblem&;

//     /// Add a given amount of a species to the equilibrium recipe
//     /// The species must be present in the chemical system.
//     /// @param name The name of the species
//     /// @param amount The amount of the species
//     /// @param unit The unit of the amount (must be convertible to either mol or kg)
//     auto addSpecies(std::string name, double amount, std::string unit) -> EquilibriumProblem&;

//     /// Add the mole amounts of the equilibrium species in a ChemicalState instance to the equilibrium recipe.
//     /// This method only extracts the mole amounts of equilibrium species in the given chemical state.
//     /// If not all species in the system is under equilibrium assumption for this equilibrium calculation,
//     /// ensure you set the partition of the chemical system before calling this method.
//     /// @note If a multiplication factor is needed, for example `2.0`, use `addState(2.0*state)`.
//     /// @param state The ChemicalState instance with the mole amounts of the species
//     /// @see setPartition
//     auto addState(const ChemicalState& state) -> EquilibriumProblem&;

//     /// Return a reference to the ChemicalSystem instance used to create this EquilibriumProblem instance
//     auto system() const -> const ChemicalSystem&;

//     /// Return the temperature for the equilibrium calculation (in unit of K)
//     auto temperature() const -> double;

//     /// Return the pressure for the equilibrium calculation (in unit of Pa)
//     auto pressure() const -> double;

//     /// Return the amounts of the elements for the equilibrium calculation (in unit of mol)
//     auto elementAmounts() const -> VectorXrConstRef;

// private:
//     struct Impl;

//     std::unique_ptr<Impl> pimpl;
// };

// /// The auxiliary class used to define equilibrium controls.
// class EquilibriumProblem::Control
// {
// public:
//     /// Construct a default EquilibriumProblem::Control object.
//     Control();

//     /// Enable temperature control during the chemical equilibrium calculation.
//     auto temperature() -> Control&;

//     /// Enable pressure control during the chemical equilibrium calculation.
//     auto pressure() -> Control&;

//     /// Enable titration of a substance during the chemical equilibrium calculation.
//     auto titrationOf(String titrant) -> Control&;

//     // Let parent class EquilibriumProblem be friend for private access.
//     friend EquilibriumProblem;

// private:
//     /// The boolean flag that indicates whether temperature is controlled or not.
//     bool controlling_temperature;

//     /// The boolean flag that indicates whether pressure is controlled or not.
//     bool controlling_pressure;

//     /// The list of substances that need to be titrated in/out.
//     Strings controlling_titrants;

//     /// The list of substance pairs that cannot be titrated simultaneously.
//     Pairs<String, String> mutually_exclusive_titrants;
// };

// /// The auxiliary class used to define equilibrium constraints.
// class EquilibriumProblem::Until
// {
// public:
//     /// Construct a default EquilibriumProblem::Until object.
//     Until();

//     /// Enforce a value for the **enthalpy** of the system at the chemical equilibrium state.
//     /// @param value The constrained enthalpy of the system
//     /// @param unit The unit of the constrained enthalpy value (must be convertible to J)
//     auto enthalpy(real value, String unit) -> Until&;

//     /// Enforce a value for the **entropy** of the system at the chemical equilibrium state.
//     /// @param value The constrained entropy of the system
//     /// @param unit The unit of the constrained entropy value (must be convertible to J/K)
//     auto entropy(real value, String unit) -> Until&;

//     /// Enforce a value for the **Gibbs energy** of the system at the chemical equilibrium state.
//     /// @param value The constrained Gibbs energy of the system
//     /// @param unit The unit of the constrained Gibbs energy value (must be convertible to J)
//     auto gibbsEnergy(real value, String unit) -> Until&;

//     /// Enforce a value for the **Helmholtz energy** of the system at the chemical equilibrium state.
//     /// @param value The constrained Helmholtz energy of the system
//     /// @param unit The unit of the constrained Helmholtz energy value (must be convertible to J)
//     auto helmholtzEnergy(real value, String unit) -> Until&;

//     /// Enforce a value for the **internal energy** of the system at the chemical equilibrium state.
//     /// @param value The constrained internal energy of the system
//     /// @param unit The unit of the constrained internal energy value (must be convertible to J)
//     auto internalEnergy(real value, String unit) -> Until&;

//     /// Enforce a value for the **volume** of the system at the chemical equilibrium state.
//     /// @param value The constrained volume of the system
//     /// @param unit The unit of the constrained volume value (must be convertible to m@sup{3})
//     auto volume(real value, String unit) -> Until&;

//     /// Enforce a value for the **chemical potential** of a species at the chemical equilibrium state.
//     /// @param species The name of the chemical species as found in the database under use
//     /// @param value The constrained chemical potential value
//     /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
//     /// @note The chemical species need not be in the chemical system; only in the database.
//     /// @warning An error will be raised if the database does not contain a species with given name.
//     auto chemicalPotential(String species, real value, String unit) -> Until&;

//     /// Enforce a value for the **chemical potential** of a species at the chemical equilibrium state.
//     /// @param species The name of the chemical species as found in the database under use
//     /// @param fn The function that computes the constrained chemical potential value at given temperature and pressure
//     /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
//     /// @note The chemical species need not be in the chemical system; only in the database.
//     /// @warning An error will be raised if the database does not contain a species with given name.
//     auto chemicalPotential(String species, Fn<real(real,real)> fn, String unit) -> Until&;

//     /// Enforce a value for the **activity** of a species at the chemical equilibrium state.
//     /// @param species The name of the chemical species as found in the database under use
//     /// @param value The constrained activity value
//     /// @note The chemical species need not be in the chemical system; only in the database.
//     /// @warning An error will be raised if the database does not contain a species with given name.
//     auto activity(String species, real value) -> Until&;

//     /// Enforce a value for the **fugacity** of a gaseous species at the chemical equilibrium state.
//     /// @param species The name of the gaseous species as found in the database under use
//     /// @param value The constrained fugacity value
//     /// @param unit The unit for the constrained fugacity value (must be convertible to Pa)
//     /// @note The gaseous species need not be in the chemical system; only in the database.
//     /// @warning An error will be raised if the database does not contain a species with given name.
//     /// @warning It is also not permitted to use a species that is not a gas.
//     auto fugacity(String species, real value, String unit) -> Until&;

//     /// Enforce a value for pH at the chemical equilibrium state.
//     /// @param value The constrained value for pH
//     auto pH(real value) -> Until&;

//     /// Enforce a value for pe at the chemical equilibrium state.
//     /// @param value The constrained value for pe
//     auto pe(real value) -> Until&;

//     /// Enforce a value for Eh at the chemical equilibrium state.
//     /// @param value The constrained value for Eh
//     /// @param unit The unit of the constrained value for Eh (must be convertible to V)
//     auto Eh(real value, String unit) -> Until&;

//     // Let parent class EquilibriumProblem be friend for private access.
//     friend EquilibriumProblem;

// private:

// };


// } // namespace Reaktoro
