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

// C++ includes
#include <string>
#include <vector>

namespace Reaktoro {
namespace kwd {

/// A type used to represent a pair (value, units).
struct ValueUnits
{
    /// Construct a default ValueUnits instance.
    ValueUnits();

    /// Construct a ValueUnits instance with given value and units.
    ValueUnits(double value, std::string units);

    /// Construct a ValueUnits instance with formatted string.
    /// The string must be in a format `value units`.
    ValueUnits(std::string str);

    /// The value
    double value;

    /// The units
    std::string units;
};

/// A type used to represent a pair (entity, value).
struct EntityValue
{
    /// Construct a default EntityValue instance.
    EntityValue();

    /// Construct a EntityValue instance with given entity, value, and units.
    EntityValue(std::string entity, double value, std::string units);

    /// Construct a EntityValue instance with formatted string.
    /// The string must be in a format `entity value units`.
    EntityValue(std::string str);

    /// The entity name
    std::string entity;

    /// The value
    double value;
};

/// A type used to represent a triplet (entity, value, units).
struct EntityValueUnits : ValueUnits
{
    /// Construct a default EntityValueUnits instance.
    EntityValueUnits();

    /// Construct a EntityValueUnits instance with given entity, value, and units.
    EntityValueUnits(std::string entity, double value, std::string units);

    /// Construct a EntityValueUnits instance with formatted string.
    /// The string must be in a format `entity value units`.
    EntityValueUnits(std::string str);

    /// The entity name
    std::string entity;
};

/// A type used to represent a triplet (value, units, entity).
struct ValueUnitsEntity : EntityValueUnits
{
    /// Construct a default ValueUnitsEntity instance.
    ValueUnitsEntity();

    /// Construct a EntityValueUnits instance with given value, units, and entity.
    ValueUnitsEntity(double value, std::string units, std::string entity);

    /// Construct a ValueUnitsEntity instance with formatted string.
    /// The string must be in a format `value units entity`.
    ValueUnitsEntity(std::string str);
};

/// A type used to represent a list of gases.
struct Gases : std::vector<std::string> {};

/// A type used to represent a list of pure minera phases.
struct Minerals: std::vector<std::string> {};

/// A type used to represent a recipe keyword.
/// A recipe keyword is used to define a list of compounds and their amounts
/// as triplets `amount units compound` (e.g., `1 kg H2O`, `1 mmol NaCl`).
struct Recipe : std::vector<ValueUnitsEntity>
{
};

/// A type used to represent keywords for equilibrium constraints.
struct EquilibriumConstraintBase : EntityValueUnits
{
    /// The titrant name used to control the constraint (if any).
    std::string titrant1;

    /// The additional titrant name used to control the constraint (if any).
    /// If titrant1 and titrant2 are non-empty, then they are mutually exclusive,
    /// which means only one will exist in positive amounts, while the other is zero.
    std::string titrant2;
};

struct pH              : EquilibriumConstraintBase {};
struct SpeciesAmount   : EquilibriumConstraintBase {};
struct SpeciesActivity : EquilibriumConstraintBase {};
struct SpeciesFugacity : EquilibriumConstraintBase {};
struct PhaseAmount     : EquilibriumConstraintBase {};
struct PhaseVolume     : EquilibriumConstraintBase {};

/// A type used to represent a chemical plot.
struct Plot
{
    /// The name of the plot file.
    std::string name;

    /// The quantity to be plot along the x-axis.
    std::string x;

    /// The quantities to be plot along the y-axis.
    std::string y;

    /// The label used for the x-axis.
    std::string xlabel;

    /// The label used for the y-axis.
    std::string ylabel;

    /// The titles of the data plot along the y-axis.
    std::string ytitles;

    /// The settings for the key.
    std::string key;
};

/// A type used to represent an equilibrium calculation.
struct EquilibriumProblem
{
    /// The name of the chemical state where this equilibrium calculation is saved.
    std::string stateid = "State";

    /// The temperature for the equilibrium calculation.
    ValueUnits temperature = {25.0, "celsius"};

    /// The pressure for the equilibrium calculation.
    ValueUnits pressure = {1.0, "bar"};

    /// The recipe definition for the equilibrium calculation.
    Recipe recipe;

    /// The pH constraints (only the last one used)
    std::vector<pH> ph;

    /// The species amount constraints
    std::vector<SpeciesAmount> species_amounts;

    /// The species activity constraints
    std::vector<SpeciesActivity> species_activities;

    /// The species fugacity constraints
    std::vector<SpeciesFugacity> species_fugacities;

    /// The phase amount constraints
    std::vector<PhaseAmount> phase_amounts;

    /// The phase volume constraints
    std::vector<PhaseVolume> phase_volumes;

    /// The list of inert species in the equilibrium calculation.
    std::vector<EntityValueUnits> inert_species;

    /// The list of inert phases in the equilibrium calculation.
    std::vector<std::string> inert_phases;
};

/// A type used to represent an equilibrium path calculation.
struct EquilibriumPath
{
    /// The name of the chemical state describing the initial condition for the equilibrium path.
    std::string initial_state;

    /// The name of the chemical state describing the final condition for the equilibrium path.
    std::string final_state;

    /// The names of the species that are hold inert in the calculation.
    std::vector<std::string> inert_species;

    /// The plots to be executed during the calculation.
    std::vector<Plot> plots;
};

/// A type used to represent a kinetic path calculation.
struct KineticPath
{
    /// The name of the chemical state where this kinetic calculation is saved.
    std::string stateid = "State";

    /// The name of the initial chemical state from where this kinetic calculation should start.
    std::string initial_condition = "State";

    /// The names of the species that are hold inert in the calculation.
    std::vector<std::string> inert_species;

    /// The names of the species that are controlled by kinetics
    std::vector<std::string> kinetic_species;

    /// The duration of the kinetic calculation
    ValueUnits duration;

    /// The plots to be executed during the calculation.
    std::vector<Plot> plots;
};

/// A type used to represent a mineral reaction.
struct MineralReaction
{
    /// The name of the mineral for this reaction.
    std::string mineral;

    /// The equation of this mineral reaction.
    std::string equation;

    /// The kinetic mechanisms of this mineral reaction.
    std::vector<std::string> mechanisms;

    /// The specific surface area of this mineral.
    ValueUnits ssa;
};

/// A type used to represent the list of compounds and their concentrations
/// in a speciation calculation.
/// The Concentrations keyword is used to define a list of compounds and
/// their concentrations as pairs `compound concentration` (e.g., `NaCl 1;
/// CaCl2 2`). Units are not needed, as they are set as a common concentration
/// units for all inputs.
struct Concentrations : std::vector<EntityValue>
{
};

/// A type used to represent an equilibrium speciation calculation.
struct SpeciationProblem
{
    /// The name of the chemical state where this equilibrium calculation is saved.
    std::string stateid = "Speciation";

    /// The temperature for the equilibrium calculation.
    ValueUnits temperature = {25.0, "celsius"};

    /// The pressure for the equilibrium calculation.
    ValueUnits pressure = {1.0, "bar"};

    /// The common units for the speciation input.
    std::string units;

    /// The concentration inputs for the speciation calculation.
    Concentrations concentrations;

    /// The pH constraints (only the last one used)
    std::vector<pH> ph;

    /// The species activity constraints
    std::vector<SpeciesActivity> species_activities;

    /// The species fugacity constraints
    std::vector<SpeciesFugacity> species_fugacities;

    /// The list of inert species in the equilibrium calculation.
    std::vector<EntityValueUnits> inert_species;

    /// The list of inert phases in the equilibrium calculation.
    std::vector<std::string> inert_phases;
};

/// A type used to represent a Phreeqc keyword.
struct PhreeqcKeyword
{
	/// The name of the chemical state after PHREEQC input script is executed.
    std::string stateid;

    /// The Phreeqc database file name.
    std::string database = "phreeqc.dat";

    /// The Phreeqc input file name.
    std::string input;

    /// The Phreeqc output file name.
    std::string output;
};

} // namespace kwd
} // namespace Reaktoro
