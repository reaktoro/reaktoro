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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;

/// The class used to define conditions to be satisfied at chemical equilibrium.
class EquilibriumConditions
{
public:
    /// Construct an EquilibriumConditions object.
    explicit EquilibriumConditions(EquilibriumSpecs const& specs);

    //=================================================================================================
    //
    // METHODS TO SPECIFY THERMODYNAMIC CONDITIONS
    //
    //=================================================================================================

    /// Specify the **temperature** of the system at chemical equilibrium.
    /// @param value The temperature of the system
    /// @param unit The unit of the temperature value (must be convertible to K)
    auto temperature(real const& value, String const& unit="K") -> void;

    /// Specify the **pressure** of the system at chemical equilibrium.
    /// @param value The pressure of the system
    /// @param unit The unit of the pressure value (must be convertible to Pa)
    auto pressure(real const& value, String const& unit="Pa") -> void;

    /// Specify the **volume** of the system at chemical equilibrium.
    /// @param value The volume of the system
    /// @param unit The unit of the volume value (must be convertible to m@sup{3})
    auto volume(real const& value, String const& unit="m3") -> void;

    /// Specify the **internal energy** of the system at chemical equilibrium.
    /// @param value The internal energy of the system
    /// @param unit The unit of the internal energy value (must be convertible to J)
    auto internalEnergy(real const& value, String const& unit="J") -> void;

    /// Specify the **enthalpy** of the system at chemical equilibrium.
    /// @param value The enthalpy of the system
    /// @param unit The unit of the enthalpy value (must be convertible to J)
    auto enthalpy(real const& value, String const& unit="J") -> void;

    /// Specify the **Gibbs energy** of the system at chemical equilibrium.
    /// @param value The Gibbs energy of the system
    /// @param unit The unit of the Gibbs energy value (must be convertible to J)
    auto gibbsEnergy(real const& value, String const& unit="J") -> void;

    /// Specify the **Helmholtz energy** of the system at chemical equilibrium.
    /// @param value The Helmholtz energy of the system
    /// @param unit The unit of the Helmholtz energy value (must be convertible to J)
    auto helmholtzEnergy(real const& value, String const& unit="J") -> void;

    /// Specify the **entropy** of the system at chemical equilibrium.
    /// @param value The entropy of the system
    /// @param unit The unit of the entropy value (must be convertible to J/K)
    auto entropy(real const& value, String const& unit="J/K") -> void;

    /// Specify the **electric charge** at chemical equilibrium.
    /// @param value The electric charge amount in the system
    /// @param unit The unit of the electric charge value (must be convertible to mol)
    auto charge(real const& value, String const& unit="mol") -> void;

    /// Specify the **amount of an element** at chemical equilibrium.
    /// @param element The name or index of the element in the system
    /// @param value The amount of an element in the system
    /// @param unit The unit of the element amount value (must be convertible to mol)
    auto elementAmount(StringOrIndex const& element, real const& value, String const& unit="mol") -> void;

    /// Specify the **amount of an element in a phase** at chemical equilibrium.
    /// @param element The name or index of the element in the system
    /// @param phase The name or index of the phase in the system
    /// @param value The amount of an element in a phase
    /// @param unit The unit of the element amount value (must be convertible to mol)
    auto elementAmountInPhase(StringOrIndex const& element, StringOrIndex const& phase, real const& value, String const& unit="mol") -> void;

    /// Specify the **mass of an element** at chemical equilibrium.
    /// @param element The name or index of the element in the system
    /// @param value The mass of an element in the system
    /// @param unit The unit of the element mass value (must be convertible to kg)
    auto elementMass(StringOrIndex const& element, real const& value, String const& unit="kg") -> void;

    /// Specify the **mass of an element in a phase** at chemical equilibrium.
    /// @param element The name or index of the element in the system
    /// @param phase The name or index of the phase in the system
    /// @param value The mass of an element in a phase
    /// @param unit The unit of the element mass value (must be convertible to kg)
    auto elementMassInPhase(StringOrIndex const& element, StringOrIndex const& phase, real const& value, String const& unit="kg") -> void;

    /// Specify the **amount of a phase** at chemical equilibrium.
    /// @param phase The name or index of the phase in the system
    /// @param value The amount of a phase in the system
    /// @param unit The unit of the phase amount value (must be convertible to mol)
    auto phaseAmount(StringOrIndex const& phase, real const& value, String const& unit="mol") -> void;

    /// Specify the **mass of a phase** at chemical equilibrium.
    /// @param phase The name or index of the phase in the system
    /// @param value The mass of a phase in the system
    /// @param unit The unit of the phase mass value (must be convertible to kg)
    auto phaseMass(StringOrIndex const& phase, real const& value, String const& unit="kg") -> void;

    /// Specify the **volume of a phase** at chemical equilibrium.
    /// @param phase The name or index of the phase in the system
    /// @param value The volume of a phase in the system
    /// @param unit The unit of the phase volume value (must be convertible to m3)
    auto phaseVolume(StringOrIndex const& phase, real const& value, String const& unit="m3") -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY SURFACE AREA CONDITIONS
    //
    //=================================================================================================

    /// Specify the **surface areas of all reacting phase interfaces** at chemical equilibrium.
    /// The surface areas that are unknown in the calculation are ignored in this method.
    /// @param values The area of all surfaces in the system (in m2)
    auto surfaceAreas(ArrayXrConstRef const& values) -> void;

    /// Specify the **surface area of a reacting phase interface** at chemical equilibrium.
    /// @param surface The name or index of the surface in the system
    /// @param value The area of a surface in the system
    /// @param unit The unit of the surface area value (must be convertible to m2)
    auto surfaceArea(StringOrIndex const& surface, real const& value, String const& unit="m2") -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY CHEMICAL POTENTIAL CONDITIONS
    //
    //=================================================================================================

    /// Specify the **chemical potential** of a substance at chemical equilibrium.
    /// @param substance The chemical formula of the substance.
    /// @param value The constrained chemical potential value.
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a chemical potential constraint for the substance.
    auto chemicalPotential(String const& substance, real const& value, String const& unit="J/mol") -> void;

    /// Specify the **ln activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained ln activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto lnActivity(String const& species, real const& value) -> void;

    /// Specify the **lg activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained lg activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto lgActivity(String const& species, real const& value) -> void;

    /// Specify the **activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto activity(String const& species, real const& value) -> void;

    /// Specify the **fugacity** of a gaseous species at chemical equilibrium.
    /// @param species The name of the gaseous species.
    /// @param value The constrained fugacity value.
    /// @param unit The unit for the constrained fugacity value (must be convertible to Pa; default is bar).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a fugacity constraint for the gas.
    auto fugacity(String const& species, real const& value, String const& unit="bar") -> void;

    /// Specify the *pH* at chemical equilibrium.
    /// @param value The constrained value for pH.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pH constraint.
    auto pH(real const& value) -> void;

    /// Specify the *pMg* at chemical equilibrium.
    /// @param value The constrained value for pMg.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pMg constraint.
    auto pMg(real const& value) -> void;

    /// Specify the *pE* at chemical equilibrium.
    /// @param value The constrained value for pE.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pE constraint.
    auto pE(real const& value) -> void;

    /// Specify the *Eh* at chemical equilibrium.
    /// @param value The constrained value for Eh.
    /// @param unit The unit of the constrained value for Eh (must be convertible to V).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an Eh constraint.
    auto Eh(real const& value, String const& unit="V") -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY LOWER AND UPPER BOUNDS FOR UNKNOWN VARIABLES
    //
    //=================================================================================================

    /// Set the lower bound for temperature during the equilibrium calculation.
    auto setLowerBoundTemperature(double value, String const& unit="K") -> void;

    /// Set the upper bound for temperature during the equilibrium calculation.
    auto setUpperBoundTemperature(double value, String const& unit="K") -> void;

    /// Set the lower bound for pressure during the equilibrium calculation.
    auto setLowerBoundPressure(double value, String const& unit="Pa") -> void;

    /// Set the upper bound for pressure during the equilibrium calculation.
    auto setUpperBoundPressure(double value, String const& unit="Pa") -> void;

    /// Set the lower bound for the amount of a titrant during the equilibrium calculation.
    auto setLowerBoundTitrant(String const& substance, double value, String const& unit="mol") -> void;

    /// Set the upper bound for the amount of a titrant during the equilibrium calculation.
    auto setUpperBoundTitrant(String const& substance, double value, String const& unit="mol") -> void;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Set the value of an input variable with given name.
    /// @param name The unique name of the input variable.
    /// @param value The new value of the input variable.
    /// @warning An error is thrown if there are no input variable with given name.
    auto set(String const& input, real const& val) -> void;

    /// Set the value of an input variable with given name.
    /// @param name The unique name of the input variable.
    /// @param value The new value of the input variable.
    /// @warning An error is thrown if there are no input variable with given name.
    auto setInputVariable(String const& name, real const& val) -> void;

    /// Set the value of an input variable with given index.
    /// @param index The index of the input variable.
    /// @param value The new value of the input variable.
    /// @warning An error is thrown if given index is out of bounds.
    auto setInputVariable(Index index, real const& val) -> void;

    /// Set the input variables with given vector of input values.
    auto setInputVariables(ArrayXrConstRef const& values) -> void;

    /// Set the values of the specified lower bounds for the *p* control variables.
    auto setLowerBoundsControlVariablesP(ArrayXdConstRef const& values) -> void;

    /// Set the values of the specified upper bounds for the *p* control variables.
    auto setUpperBoundsControlVariablesP(ArrayXdConstRef const& values) -> void;

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> ChemicalSystem const&;

    /// Return the names of the input variables associated with the equilibrium conditions.
    auto inputNames() const -> Strings const&;

    /// Return the values of the input variables associated with the equilibrium conditions.
    auto inputValues() const -> ArrayXrConstRef;

    /// Return the values of the input variables associated with the equilibrium conditions.
    /// @param name The unique name of the *w* input variable
    auto inputValue(String const& name) const -> real const&;

    /// Return the specified lower bounds for the *p* control variables.
    auto lowerBoundsControlVariablesP() const -> ArrayXdConstRef;

    /// Return the specified upper bounds for the *p* control variables.
    auto upperBoundsControlVariablesP() const -> ArrayXdConstRef;

private:
    const ChemicalSystem msystem;         ///< The chemical system associated with these equilibrium conditions.
    const Strings wvars;                  ///< The names of the *w* input variables in the equilibrium problem.
    const Strings pvars;                  ///< The names of the *p* control variables variables in the equilibrium problem.
    const Index itemperature_p;           ///< The index of the temperature variable among the *p* control variables.
    const Index ipressure_p;              ///< The index of the pressure variable among the *p* control variables.
    const Indices isurface_areas_w;       ///< The indices of the surface area variables among the *w* input variables.
    const Indices isurface_areas_p;       ///< The indices of the surface area variables among the *p* control variables.
    const Indices isurface_areas_known;   ///< The indices of the surfaces areas that are known.
    const Indices isurface_areas_unknown; ///< The indices of the surfaces areas that are unknown.
    ArrayXr w;                           ///< The values of the *w* input variables.
    ArrayXd plower;                      ///< The lower bounds for the *p* control variables.
    ArrayXd pupper;                      ///< The upper bounds for the *p* control variables.
};

} // namespace Reaktoro
