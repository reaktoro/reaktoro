// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class EquilibriumSpecs;
class Param;

/// The sensitivity derivatives of a chemical equilibrium state.
/// This class stores the sensitivity derivatives of a chemical equilibrium
/// state. These are partial derivatives of the species amounts with respect to
/// input variables, such as temperature, pressure, amounts of the elements.
/// If the chemical equilibrium state is computed with other given conditions,
/// for example, given volume and internal energy, derivatives with respect to
/// these input conditions will be available. These sensitivity derivatives are
/// important for implicit numerical methods, since they enable faster
/// convergence rates.
class EquilibriumSensitivity
{
public:
    /// Construct a default EquilibriumSensitivity object.
    EquilibriumSensitivity();

    /// Construct an EquilibriumSensitivity object with given equilibrium problem specifications.
    EquilibriumSensitivity(const EquilibriumSpecs& specs);

    /// Construct a copy of an EquilibriumSensitivity object.
    EquilibriumSensitivity(const EquilibriumSensitivity& other);

    /// Destroy this EquilibriumSensitivity object.
    ~EquilibriumSensitivity();

    /// Assign a copy of an EquilibriumSensitivity object to this.
    auto operator=(EquilibriumSensitivity other) -> EquilibriumSensitivity&;

    /// Initialize this EquilibriumSensitivity object with given equilibrium problem specifications.
    auto initialize(const EquilibriumSpecs& specs) -> void;

    //======================================================================
    // DERIVATIVES OF SPECIES AMOUNTS WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the species amounts *n* with respect to an input variable in *w*.
    /// @param wid The identifier of the input variable in *w* (e.g., "T", "P", "pH", it depends on what is input).
    auto dndw(const String& wid) const -> VectorXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to an input variable in *w*.
    /// @param param The input variable in *w* as a Param object.
    auto dndw(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to the input variables *w*.
    auto dndw() const -> MatrixXdConstRef;

    /// Set the derivatives of the species amounts *n* with respect to the input variables *w*.
    auto dndw(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF p-CONTROL VARIABLES WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the *p* control variables with respect to an input variable in *w*.
    /// @param wid The identifier of the input variable in *w* (e.g., "T", "P", "pH", it depends on what is input).
    auto dpdw(const String& wid) const -> VectorXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to an input variable in *w*.
    /// @param param The input variable in *w* as a Param object.
    auto dpdw(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to the input variables *w*.
    auto dpdw() const -> MatrixXdConstRef;

    /// Set the derivatives of the *p* control variables with respect to the input variables *w*.
    auto dpdw(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF q-CONTROL VARIABLES WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the *q* control variables with respect to an input variable in *w*.
    /// @param wid The identifier of the input variable in *w* (e.g., "T", "P", "pH", it depends on what is input).
    auto dqdw(const String& wid) const -> VectorXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to an input variable in *w*.
    /// @param param The input variable in *w* as a Param object.
    auto dqdw(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to the input variables *w*.
    auto dqdw() const -> MatrixXdConstRef;

    /// Set the derivatives of the *q* control variables with respect to the input variables *w*.
    auto dqdw(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF SPECIES AMOUNTS WITH RESPECT TO COMPONENT AMOUNTS
    //======================================================================

    /// Return the derivatives of the species amounts *n* with respect to component amounts *b*.
    auto dndb() const -> MatrixXdConstRef;

    /// Set the derivatives of the species amounts *n* with respect to component amounts *b*.
    auto dndb(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF p-CONTROL VARIABLES WITH RESPECT TO COMPONENT AMOUNTS
    //======================================================================

    /// Return the derivatives of the control variables *p* with respect to component amounts *b*.
    auto dpdb() const -> MatrixXdConstRef;

    /// Set the derivatives of the control variables *p* with respect to component amounts *b*.
    auto dpdb(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF q-CONTROL VARIABLES WITH RESPECT TO COMPONENT AMOUNTS
    //======================================================================

    /// Return the derivatives of the control variables *q* with respect to component amounts *b*.
    auto dqdb() const -> MatrixXdConstRef;

    /// Set the derivatives of the control variables *q* with respect to component amounts *b*.
    auto dqdb(MatrixXdConstRef data) -> void;

    //======================================================================
    // TOTAL DERIVATIVES OF CHEMICAL PROPERTIES
    //======================================================================

    /// Return the total derivatives of the chemical properties *u* with respect to input variables *w*.
    auto dudw() const -> MatrixXdConstRef;

    /// Return the total derivatives of the chemical properties *u* with respect to component amounts *b*.
    auto dudb() const -> MatrixXdConstRef;

    /// Set the total derivatives of the chemical properties *u* with respect to input variables *w*.
    auto dudw(MatrixXdConstRef data) -> void;

    /// Set the total derivatives of the chemical properties *u* with respect to component amounts *b*.
    auto dudb(MatrixXdConstRef data) -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
