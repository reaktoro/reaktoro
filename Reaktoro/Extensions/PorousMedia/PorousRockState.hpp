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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalSystem;
class ChemicalState;

// TODO: Implement PorousRockState class in the PorousMedia extension.

/// The chemical state of a porous rock system filled with fluids.
/// @ingroup PorousMediaExtension
class PorousRockState : public ChemicalState
{
public:
    /// Construct a PorousRockState instance with standard conditions.
    /// This constructor creates an instance of PorousRockState with temperature
    /// 25 °C, pressure 1 bar, and zero mole amounts for the species.
    explicit PorousRockState(const ChemicalSystem& system);

    /// Construct a copy of a PorousRockState instance.
    PorousRockState(const PorousRockState& other);

    /// Destroy this PorousRockState instance.
    virtual ~PorousRockState();

    /// Assign a PorousRockState instance to this instance.
    auto operator=(PorousRockState other) -> PorousRockState&;

    /// Set the porosity of the rock system.
    auto setPorosity(real value) -> void;

    /// Set the volume fraction of the fluid phase with given index among all other fluid phases.
    auto setVolumeFractionAmongFluids(Index iphase, real value) -> void;

    /// Set the volume fraction of the fluid phase with given named among all other fluid phases.
    auto setVolumeFractionAmongFluids(String phasename, real value) -> void;

    /// Set the volume fraction of the solid phase with given index among all other solid phases.
    auto setVolumeFractionAmongSolids(Index iphase, real value) -> void;

    /// Set the volume fraction of the solid phase with given named among all other solid phases.
    auto setVolumeFractionAmongSolids(String phasename, real value) -> void;

    /// Perform an equilibrium calculation of the porous rock system.
    auto equilibrate() -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output a PorousRockState instance.
auto operator<<(std::ostream& out, const PorousRockState& state) -> std::ostream&;

} // namespace Reaktoro
