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
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>

namespace Reaktoro {

/// The sensitivity derivatives of a chemical equilibrium state.
/// This class stores the sensitivity derivatives of a chemical equilibrium
/// state. These are partial derivatives of the species amounts with respect to
/// input variables, such as temperature, pressure, amounts of the elements.
/// If the chemical equilibrium state is computed with other given conditions,
/// for example, given volume and internal energy, derivatives with respect to
/// these input conditions will be available. These sensitivity derivatives are
/// important for implicit numerical methods, since they enable faster
/// convergence rates.
class KineticsSensitivity : public EquilibriumSensitivity
{
public:
    /// Construct a default KineticsSensitivity object.
    KineticsSensitivity();

    /// Construct a KineticsSensitivity object with given equilibrium problem specifications.
    explicit KineticsSensitivity(EquilibriumSpecs const& specs);

    /// Construct a KineticsSensitivity object from an EquilibriumSensitivity object.
    KineticsSensitivity(EquilibriumSensitivity const& other);

    /// Initialize this KineticsSensitivity object with given equilibrium problem specifications.
    auto initialize(EquilibriumSpecs const& specs) -> void;

    //============================================================================================
    // SENSITIVITY DERIVATIVES OF SPECIES AMOUNTS WITH RESPECT TO SURFACE AREAS AND TIME STEP
    //============================================================================================

    /// Return the derivatives of the species amounts *n* with respect to a surface area in vector *S*.
    /// @param sid The identifier of the reacting interphase surface (e.g., "AqueousPhase:GaseousPhase", "Calcite").
    auto dnds(String const& name) const -> VectorXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to the surface areas in vector *S*.
    auto dnds() const -> MatrixXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to time step *Δt*.
    auto dnddt() const -> VectorXdConstRef;

    //============================================================================================
    // SENSITIVITY DERIVATIVES OF p-CONTROL VARIABLES WITH RESPECT TO SURFACE AREAS AND TIME STEP
    //============================================================================================

    /// Return the derivatives of the *p* control variables with respect to a surface area in vector *S*.
    /// @param sid The identifier of the reacting interphase surface (e.g., "AqueousPhase:GaseousPhase", "Calcite").
    auto dpds(String const& name) const -> VectorXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to the surface areas in vector *S*.
    auto dpds() const -> MatrixXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to time step *Δt*.
    auto dpddt() const -> VectorXdConstRef;

    //============================================================================================
    // SENSITIVITY DERIVATIVES OF q-CONTROL VARIABLES WITH RESPECT TO SURFACE AREAS AND TIME STEP
    //============================================================================================

    /// Return the derivatives of the *q* control variables with respect to a surface area in vector *S*.
    /// @param sid The identifier of the reacting interphase surface (e.g., "AqueousPhase:GaseousPhase", "Calcite").
    auto dqds(String const& name) const -> VectorXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to the input variables *w*.
    auto dqds() const -> MatrixXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to time step *Δt*.
    auto dqddt() const -> VectorXdConstRef;

    //============================================================================================
    // SENSITIVITY DERIVATIVES OF CHEMICAL PROPERTIES WITH RESPECT TO SURFACE AREAS AND TIME STEP
    //============================================================================================

    /// Return the sensitivity derivatives of the chemical properties *u* with respect to the surface areas in vector *S*.
    auto duds() const -> MatrixXdConstRef;

    /// Return the sensitivity derivatives of the chemical properties *u* with respect to time step *Δt*.
    auto duddt() const -> VectorXdConstRef;
};

} // namespace Reaktoro
