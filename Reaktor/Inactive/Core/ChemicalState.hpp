/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <memory>
#include <string>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/VectorResult.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalSystem;
class Partitioning;

/**
 * Provides a computational representation of the state of a multiphase chemical system
 *
 * The chemical state of a multiphase system is defined by its temperature @f$(T)@f$,
 * pressure @f$(P)@f$, and molar composition @f$(\mathbf{n})@f$.
 *
 * **Usage**
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * // Create a ChemicalState instance, where system is a ChemicalSystem instance
 * ChemicalState state(system);
 *
 * // Set the temperature and pressure states
 * state.setTemperature(60.0, "celsius");
 * state.setPressure(  180.0, "bar");
 *
 * // Set the amount of some species
 * state.set( "H2O(l)",  1.0, "kg");
 * state.set(    "Na+",  1.0, "mol");
 * state.set(    "Cl-",  1.0, "mol");
 * state.set("CO2(aq)",  0.5, "mol");
 * state.set("Calcite", 10.0, "g");
 *
 * // Output the chemical state instance
 * std::cout << state << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @see ChemicalSystem
 *
 * @ingroup Core
 */
class ChemicalState
{
public:
	/**
	 * Constructs a ChemicalState instance from a ChemicalSystem instance
	 * @param system The chemical system instance
	 */
	explicit ChemicalState(const ChemicalSystem& system);

	/**
	 * Constructs a copy of a ChemicalState instance
	 */
	ChemicalState(const ChemicalState& other);

	/**
	 * Destroys the instance
	 */
	virtual ~ChemicalState();

	/**
	 * Assings a ChemicalState instance to this instance
	 */
    auto operator=(ChemicalState other) -> ChemicalState&;

    /**
     * Gets the temperature of the chemical state (in units of K)
     */
    auto temperature() const -> double;

    /**
     * Gets the pressure of the chemical state (in units of Pa)
     */
    auto pressure() const -> double;

    /**
     * Gets the molar composition of the chemical state (in units of mol)
     */
    auto composition() const -> const Vector&;

    /**
     * Gets the chemical system instance
     */
    auto system() const -> const ChemicalSystem&;

    /**
     * Sets the temperature of the chemical state
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * state.setTemperature(50.0 * unit(degC));
     * state.setTemperature(323.15 * unit(K));
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param value The temperature value with unit types
     */
    auto setTemperature(units::Temperature value) -> void;

    /**
     * Sets the temperature of the chemical state
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * state.setTemperature(50.0, "degC");
     * state.setTemperature(323.15, "K");
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param value The temperature value
     * @param unit The string representing the temperature unit
     */
    auto setTemperature(double value, const std::string& unit) -> void;

	/**
     * Sets the pressure of the chemical state
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * state.setPressure(150.0 * unit(bar));
     * state.setPressure(10.0 * unit(MPa));
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param value The pressure value with unit types
     */
    auto setPressure(units::Pressure value) -> void;

    /**
     * Sets the pressure of the chemical state
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * state.setPressure(150.0, "bar");
     * state.setPressure(10.0, "MPa");
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param value The pressure value
     * @param unit The string representing the pressure unit
     */
    auto setPressure(double value, const std::string& unit) -> void;

    /**
     * Sets the molar composition of the species
     * @param value The number of moles of each species
     */
    auto setComposition(double value) -> void;

    /**
     * Sets the molar composition of the species
     * @param n The vector of number of moles of the species
     */
    auto setComposition(const Vector& n) -> void;

    /**
     * Sets the molar composition of a set of species
     * @param indices The indices of the set of species
     * @param values The molar amounts of the species in the set
     */
    auto setComposition(const Indices& indices, const Vector& values) -> void;

    /**
     * Sets the molar composition of the equilibrium species
     * @param partitioning The partitioning instance
     * @param ne The molar composition of the equilibrium species
     *
     * @see Partitioning
     */
    auto setEquilibriumComposition(const Partitioning& partitioning, const Vector& ne) -> void;

    /**
     * Sets the molar composition of the kinetic species
     * @param partitioning The partitioning instance
     * @param nk The molar composition of the kinetic species
     *
     * @see Partitioning
     */
    auto setKineticComposition(const Partitioning& partitioning, const Vector& nk) -> void;

    /**
     * Sets the molar composition of the inert species
     * @param partitioning The partitioning instance
     * @param ni The molar composition of the inert species
     *
     * @see Partitioning
     */
    auto setInertComposition(const Partitioning& partitioning, const Vector& ni) -> void;

    /**
     * Sets the amount of a species
     * @param ispecies The index of the species
     * @param value The amount of the species with units (must be convertible to mol)
     */
    auto setSpeciesAmount(Index ispecies, units::Amount value) -> void;

    /**
     * Sets the amount of a species
     * @param species The name of the species
     * @param value The amount of the species with units (must be convertible to mol)
     */
    auto setSpeciesAmount(const std::string& species, units::Amount value) -> void;

    /**
     * Sets the mass of a species
     * @param ispecies The index of the species
     * @param value The mass of the species with units (must be convertible to gram)
     */
    auto setSpeciesMass(Index ispecies, units::Mass value) -> void;

    /**
     * Sets the mass of a species
     * @param species The name of the species
     * @param value The mass of the species with units (must be convertible to gram)
     */
    auto setSpeciesMass(const std::string& species, units::Mass value) -> void;

    /**
     * Sets the composition of the mineral phases
     *
     * This method provides a convenient way to setup the molar composition of the
     * minerals from their *volume*, *porosity* and the *volume percent* of each mineral.
     *
     * **Usage**
     *
     * Assume a carbonate rock whose composition is given in the table below.
     *
     * Mineral | Volume (%)
     * :-------|:---------:
     * Calcite | 93.3
     * Dolomite| 5.2
     * Quartz  | 1.5
     *
     * Assuming that our chemical system contains 100 cm<sup>3</sup> of this rock,
     * with a porosity of 0.30, we can specify this as follows:
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * state.setMinerals("93.3:Calcite 5.2:Dolomite 1.5:Quartz", 100.0, "cm3", 0.30);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * The above call will then set the molar abundance of each mineral.
     *
     * @param composition The string containing the composition of the minerals
     * @param volume The total volume of the rock
     * @param unit The unit of the volume
     * @param porosity The porosity of the rock (0 < @f$\phi@f$ < 1)
     */
    auto setMinerals(const std::string& composition, double volume, const std::string& unit, double porosity = 0.0) -> void;

    /**
     * Sets the amount or mass of a species
     * @param ispecies The index of the species
     * @param value The amount or mass of the species (must be convertible to either mol or gram)
     */
    template<typename Unit>
    auto set(Index ispecies, units::Constant<Unit> value) -> void;

    /**
     * Sets the amount or mass of a species
     * @param species The name of the species
     * @param value The amount or mass of the species (must be convertible to either mol or gram)
     */
    template<typename Unit>
    auto set(const std::string& species, units::Constant<Unit> value) -> void;

    /**
     * Sets the amount or mass of a species
     * @param species The name of the species
     * @param value The amount or mass of the species
     * @param unit The units of the value(must be convertible to either mol or gram)
     */
    auto set(const std::string& species, double value, const std::string& unit) -> void;

    /**
     * Gets the molar composition of a phase
     * @param iphase The index of the phase
     * @return The number of moles of each species in the phase
     */
    auto phaseComposition(Index iphase) const -> Vector;

    /**
     * Gets the molar composition of the equilibrium species
     * @param partitioning The partitioning instance
     * @return The number of moles of the equilibrium species
     * @see Partitioning
     */
    auto equilibriumComposition(const Partitioning& partitioning) const -> Vector;

    /**
     * Gets the molar composition of the kinetic species
     * @param partitioning The partitioning instance
     * @return The number of moles of the kinetic species
     * @see Partitioning
     */
    auto kineticComposition(const Partitioning& partitioning) const -> Vector;

    /**
     * Gets the molar composition of the inert species
     * @param partitioning The partitioning instance
     * @return The number of moles of the inert species
     * @see Partitioning
     */
    auto inertComposition(const Partitioning& partitioning) const -> Vector;

    /**
     * Gets the amount of a species
     * @param ispecies The index of the species
     * @return The amount of the species with units
     */
    auto amount(Index ispecies) const -> units::Amount;

    /**
     * Gets the amount of a species
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * std::cout << "The number of moles of H2O(l) is ";
     * std::cout << state.amount("H2O(l)").in(unit(mol)) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param species The name of the species
     * @return The amount of the species with units
     */
    auto amount(const std::string& species) const -> units::Amount;

    /**
     * Gets the amount of an element
     * @param idx_element The index of the element
     * @return The amount of the element with units
     */
    auto amountElement(const Index& idx_element) const -> units::Amount;

    /**
     * Gets the amount of an element
     * @param element The name of the element
     * @return The amount of the element with units
     */
    auto amountElement(const std::string& element) const -> units::Amount;

    /**
     * Gets the amount of an element in a specific phase
     * @param element The name of the element
     * @param phase The name of the phase
     * @return The amount of the element in the given phase with units
     */
    auto amountElement(const std::string& element, const std::string& phase) const -> units::Amount;

    /**
     * Gets the total amount of moles in a phase
     * @param idx_phase The index of the phase
     */
    auto amountPhase(const Index& idx_phase) const -> units::Amount;

    /**
     * Gets the total amount of moles in a phase
     * @param phase The name of the phase
     */
    auto amountPhase(const std::string& phase) const -> units::Amount;

    /**
     * Gets the mass of a species
     * @param ispecies The index of the species
     * @return The mass of the species with units
     */
    auto mass(Index ispecies) const -> units::Mass;

    /**
     * Gets the mass of a species
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * std::cout << "The mass of H2O(l) is ";
     * std::cout << state.mass("H2O(l)").in(unit(kg)) << "kg" << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param species The name of the species
     * @return The mass of the species with units
     */
    auto mass(const std::string& species) const -> units::Mass;

    /**
     * Gets the mass of an element
     * @param element The name of the element
     * @return The mass of the element with units
     */
    auto massElement(const std::string& element) const -> units::Mass;

    /**
     * Gets the mass of an element in a specific phase
     * @param element The name of the element
     * @param phase The name of the phase
     * @return The mass of the element in the given phase with units
     */
    auto massElement(const std::string& element, const std::string& phase) const -> units::Mass;

    /**
     * Calculates the molality of a species
     * @param species The name of the species
     */
    auto molality(const std::string& species) const -> units::Molality;

    /**
     * Calculates the molality of an element
     * @param element The name of the element
     */
    auto molalityElement(const std::string& element) const -> units::Molality;

    /**
     * Calculates the molar fractions of the species
     */
    auto molarFractions() const -> Vector;

    /**
     * Calculates the concentrations of the species
     */
    auto concentrations() const -> Vector;

    /**
     * Calculates the activities and their molar derivatives
     */
    auto activities() const -> VectorResult;

    /**
     * Gets the activity of a species
     * @param a The activities of all species
     * @param ispecies The index of the species
     * @see activities
     */
    auto activity(const VectorResult& a, Index ispecies) const -> double;

    /**
     * Gets the activity of a species
     *
     * **Usage**
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * auto a = state.activities();
     * std::cout << "The activity of H2O(l) = " << state.activity(a, "H2O(l)") << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param a The activities of the species and their molar derivatives
     * @param ispecies The index of the species
     *
     * @see activities
     */
    auto activity(const VectorResult& a, const std::string& species) const -> double;

    /**
     * Calculates the acidity of the aqueous phase
     * @param a The activities of the species and their molar derivatives
     * @return The pH of the aqueous phase
     */
    auto acidity(const VectorResult& a) const -> double;

    /**
     * Calculates the acidity of the aqueous phase
     */
    auto acidity() const -> double;

private:
    template<typename Unit>
    auto set(Index ispecies, units::Constant<Unit> value, units::base::Mol) -> void;

    template<typename Unit>
    auto set(Index ispecies, units::Constant<Unit> value, units::base::Gram) -> void;

    template<typename Unit>
    auto set(const std::string& species, units::Constant<Unit> value, units::base::Mol) -> void;

    template<typename Unit>
    auto set(const std::string& species, units::Constant<Unit> value, units::base::Gram) -> void;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

/**
 * Outputs a @ref ChemicalState instance
 */
auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&;

template<typename Unit>
auto ChemicalState::set(Index ispecies, units::Constant<Unit> value) -> void
{
    constexpr bool value1 = units::Convertible<Unit, units::base::Mol>::value;
    constexpr bool value2 = units::Convertible<Unit, units::base::Gram>::value;
    static_assert(value1 or value2, "*** Error *** the provided unit is not convertible to units of amount or mass.");
    using unittype = typename std::conditional<value1, units::base::Mol, units::base::Gram>::type;
    set(ispecies, value, unittype());
}

template<typename Unit>
auto ChemicalState::set(const std::string& species, units::Constant<Unit> value) -> void
{
    constexpr bool value1 = units::Convertible<Unit, units::base::Mol>::value;
    constexpr bool value2 = units::Convertible<Unit, units::base::Gram>::value;
    static_assert(value1 or value2, "*** Error *** the provided unit is not convertible to units of amount or mass.");
    using unittype = typename std::conditional<value1, units::base::Mol, units::base::Gram>::type;
    set(species, value, unittype());
}

template<typename Unit>
auto ChemicalState::set(Index ispecies, units::Constant<Unit> value, units::base::Mol) -> void
{
    setSpeciesAmount(ispecies, value);
}

template<typename Unit>
auto ChemicalState::set(Index ispecies, units::Constant<Unit> value, units::base::Gram) -> void
{
    setSpeciesMass(ispecies, value);
}

template<typename Unit>
auto ChemicalState::set(const std::string& species, units::Constant<Unit> value, units::base::Mol) -> void
{
    setSpeciesAmount(species, value);
}

template<typename Unit>
auto ChemicalState::set(const std::string& species, units::Constant<Unit> value, units::base::Gram) -> void
{
    setSpeciesMass(species, value);
}

} // namespace Reaktor
