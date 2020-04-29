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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/PhaseChemicalProps.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalSystem;

// TODO: Implement AqueousProps class in the Geochemistry extension.

/// The chemical properties of an aqueous phase.
/// @ingroup GeochemistryExtension
class AqueousProps : public PhaseChemicalProps
{
public:
    /// Construct a AqueousProps instance with standard conditions.
    /// This constructor creates an instance of AqueousProps with temperature
    /// 25 °C, pressure 1 bar, zero mole amounts for the aqueous solutes, and 1
    /// kg of solvent water.
    explicit AqueousProps(const Phase& phase);

    /// Construct a copy of a AqueousProps instance.
    AqueousProps(const AqueousProps& other);

    /// Destroy this AqueousProps instance.
    virtual ~AqueousProps();

    /// Assign a AqueousProps instance to this instance.
    auto operator=(AqueousProps other) -> AqueousProps&;

    /// Return the ionic strength of the aqueous phase (in molal).
    auto ionicStrength() const -> real;

    /// Return the pH of the aqueous phase.
    auto pH() const -> real;

    /// Return the pE of the aqueous phase.
    auto pE() const -> real;

    /// Return the reduction potential of the aqueous phase (in V).
    auto Eh() const -> real;

    /// Return the total alkalinity of the aqueous phase (in eq/L).
    /// The total alkalinity (Alk) of the aqueous phase is by default
    /// calculated as the *acid neutralizing capacity* (ANC) of the solution
    /// using the formula:
    /// @f[
    /// \mathrm{Alk=[Na^{+}]+[K^{+}]+2[Ca^{2+}]+2[Mg^{2+}]-[Cl^{-}]-2[SO_{4}^{2-}]},
    /// @f] where @f$[\mathrm{species}]@f$ is the free molar concentration
    /// (mol/L) of the species in the solution. This formula is simpler,
    /// derived from the charge balance condition, and equivalent to the
    /// standard formula of total alkalinity.
    auto alkalinity() const -> real;

    /// Return the underlying Phase object for the aqueous phase.
    auto phase() const -> const Phase&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output a AqueousProps instance.
auto operator<<(std::ostream& out, const AqueousProps& state) -> std::ostream&;

} // namespace Reaktoro
