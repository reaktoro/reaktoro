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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// A class to represent a formation reaction of a chemical species.
/// @ingroup Core
class FormationReaction
{
public:
    /// Construct a default FormationReaction object.
    FormationReaction();

    /// Return a deep copy of this FormationReaction object.
    auto clone() const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new name of the product species in the formation reaction.
    auto withProduct(String product) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reactant species in the formation reaction.
    auto withReactants(Pairs<Species, double> reactants) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new equilibrium constant value (log base 10).
    auto withEquilibriumConstant(real value) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new equilibrium constant function (log base 10).
    auto withEquilibriumConstantFn(const Fn<real(real,real)>& fn) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction enthalpy value (in J/mol).
    auto withEnthalpyChange(real value) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction enthalpy function (in J/mol).
    auto withEnthalpyChangeFn(const Fn<real(real,real)>& fn) const -> FormationReaction;

    /// Return the name of the product species in the formation reaction.
    auto product() const -> String;

    /// Return the reactant species in the formation reaction.
    auto reactants() const -> const Pairs<Species, double>&;

    /// Return the equilibrium constant function of the formation reaction (log base 10).
    auto equilibriumConstantFn() const -> const Fn<real(real,real)>&;

    /// Return the reaction enthalpy function of the formation reaction.
    auto enthalpyChangeFn() const -> const Fn<real(real,real)>&;

    /// Return the standard Gibbs energy function of the product species in the formation reaction.
    auto standardGibbsEnergyFn() const -> Fn<real(real,real)>;

    /// Return the standard enthalpy function of the product species in the formation reaction.
    auto standardEnthalpyFn() const -> Fn<real(real,real)>;

    /// Return the stoichiometric coefficient of a reactant with given name in the formation reaction.
    auto stoichiometry(String reactant) const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
