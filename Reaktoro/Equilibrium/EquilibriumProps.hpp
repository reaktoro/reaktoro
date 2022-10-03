// Reaktoro is a unified framework for modeling chemically reactive phases.
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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class EquilibriumSpecs;

/// The class that computes chemical properties of a system during equilibrium calculations.
/// This class exists to allow derivatives of these chemical properties to be
/// collected during their computations.
class EquilibriumProps
{
public:
    /// Construct a EquilibriumProps object with given equilibrium specifications.
    explicit EquilibriumProps(const EquilibriumSpecs& specs);

    /// Construct a copy of a EquilibriumProps object.
    EquilibriumProps(const EquilibriumProps& other);

    /// Destroy this EquilibriumProps object.
    ~EquilibriumProps();

    /// Assign a EquilibriumProps object to this.
    auto operator=(EquilibriumProps other) -> EquilibriumProps&;

    /// Update the chemical properties of the chemical system.
    /// @param n The amounts of the species.
    /// @param p The values of the *p* control variables (e.g., T, P, n[H+] in case U, V and pH are given).
    /// @param w The input variables *w* in the chemical equilibrium problem (e.g., U, V, pH).
    /// @param useIdealModel If true, ideal thermodynamic models are used for the phases.
    auto update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel = false) -> void;

    /// Update the chemical properties of the chemical system. This is an
    /// special update method compared to ChemicalProps::update. Here, the
    /// variables in `n`, `p`, and `w` are inspected for seed state. This is
    /// seed in the sense of automatic differentiation. When one of this
    /// variables are detected to be seeded, this implies that automatic
    /// differentiation is computing derivatives of the chemical properties
    /// with respect to this seeded variable. The computed derivatives are then
    /// stored in a matrix. Access to these derivatives can be obtained with
    /// methods @ref dudn, @ref dudp, and @ref dudw.
    /// @param n The amounts of the species.
    /// @param p The values of the *p* control variables (e.g., T, P, n[H+] in case U, V and pH are given).
    /// @param w The input variables *w* in the chemical equilibrium problem (e.g., U, V, pH).
    /// @param useIdealModel If true, ideal thermodynamic models are used for the phases.
    /// @param inpw The index of the variable in (n, p, w) currently seeded for autodiff computation.
    auto update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel, long inpw) -> void;

    /// Enable recording of derivatives of the chemical properties with respect
    /// to *(n, p, w)* to contruct its full Jacobian matrix.
    /// Consider a series of forward automatic differentiation passes to
    /// compute the partial derivatives of the chemical properties with respect
    /// to the variables *(n, p, w)*. Use this method before these operations
    /// so that these derivatives are recorded in this EquilibriumProps object.
    /// At the end of all forward passes, the full Jacobian of the chemical
    /// properties will have been constructed.
    /// @note Call @ref assembleFullJacobianEnd after these forward passes
    /// have ended to eliminates the minor overhead of recording derivatives.
    /// @note Access these derivatives with methods @ref dudn, @ref dudp, and @ref dudw.
    auto assembleFullJacobianBegin() -> void;

    /// Disable recording of derivatives of the chemical properties with
    /// respect to *(n, p, w)* to indicate the end of the full Jacobian matrix
    /// construction.
    auto assembleFullJacobianEnd() -> void;

    /// Return the underlying chemical state of the system and its updated properties.
    auto chemicalState() const -> const ChemicalState&;

    /// Return the underlying chemical properties of the system.
    auto chemicalProps() const -> const ChemicalProps&;

    /// Return the partial derivatives of the serialized chemical properties
    /// *u* with respect to species amounts *n*.
    auto dudn() const -> MatrixXdConstRef;

    /// Return the partial derivatives of the serialized chemical properties
    /// *u* with respect to control variables *p*.
    auto dudp() const -> MatrixXdConstRef;

    /// Return the partial derivatives of the serialized chemical properties
    /// *u* with respect to input variables *w*.
    auto dudw() const -> MatrixXdConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
