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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class EquilibriumSpecs;
class Param;
class Params;

/// The sensitivity derivatives of a chemical equilibrium state.
/// This class stores the sensitivity derivatives of a chemical equilibrium
/// state. These are partial derivatives of the species amounts with respect to
/// input parameters, such as temperature, pressure, amounts of the elements.
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

    /// Construct a default EquilibriumSensitivity object.
    EquilibriumSensitivity(const EquilibriumSpecs& specs);

    //======================================================================
    // DERIVATIVES OF SPECIES AMOUNTS WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the species amounts *n* with respect to an input parameter in *c*.
    /// @param cid The identifier of the input parameter in *c* (e.g., "T", "P", "pH", it depends on what is input).
    auto dndc(const String& cid) const -> VectorXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to an input parameter in *c*.
    /// @param param The input parameter in *c* as a Param object.
    auto dndc(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the species amounts *n* with respect to the input parameters *c*.
    auto dndc() const -> MatrixXdConstRef;

    /// Set the derivatives of the species amounts *n* with respect to the input parameters *c*.
    auto dndc(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF p-CONTROL VARIABLES WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the *p* control variables with respect to an input parameter in *c*.
    /// @param cid The identifier of the input parameter in *c* (e.g., "T", "P", "pH", it depends on what is input).
    auto dpdc(const String& cid) const -> VectorXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to an input parameter in *c*.
    /// @param param The input parameter in *c* as a Param object.
    auto dpdc(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the *p* control variables with respect to the input parameters *c*.
    auto dpdc() const -> MatrixXdConstRef;

    /// Set the derivatives of the *p* control variables with respect to the input parameters *c*.
    auto dpdc(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF q-CONTROL VARIABLES WITH RESPECT TO INPUT PARAMETERS
    //======================================================================

    /// Return the derivatives of the *q* control variables with respect to an input parameter in *c*.
    /// @param cid The identifier of the input parameter in *c* (e.g., "T", "P", "pH", it depends on what is input).
    auto dqdc(const String& cid) const -> VectorXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to an input parameter in *c*.
    /// @param param The input parameter in *c* as a Param object.
    auto dqdc(const Param& param) const -> VectorXdConstRef;

    /// Return the derivatives of the *q* control variables with respect to the input parameters *c*.
    auto dqdc() const -> MatrixXdConstRef;

    /// Set the derivatives of the *q* control variables with respect to the input parameters *c*.
    auto dqdc(MatrixXdConstRef data) -> void;

    //======================================================================
    // DERIVATIVES OF SPECIES AMOUNTS WITH RESPECT TO COMPONENT AMOUNTS
    //======================================================================

    /// Return the derivatives of the species amounts *n* with respect to component amounts *b*.
    auto dndb() const -> MatrixXdConstRef;

    /// Set the derivatives of the species amounts *n* with respect to component amounts *b*.
    auto dndb(MatrixXdConstRef data) const -> MatrixXdConstRef;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;


};

} // namespace Reaktoro
