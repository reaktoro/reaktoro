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
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

// Forward declarations
class MineralMixture;
class MineralSpecies;

/// Class that defines an mineral phase
class MineralPhase : public Phase
{
public:
    /// Construct a default MineralPhase instance.
    MineralPhase();

    /// Construct a MineralPhase instance with given mineral mixture.
    explicit MineralPhase(const MineralMixture& mixture);

    /// Construct a MineralPhase instance with given species.
    explicit MineralPhase(const MineralSpecies& species);

    /// Set the chemical model of the phase with the ideal solution model.
    auto setChemicalModelIdeal() -> MineralPhase&;

    /// Set the chemical model of the phase with the Redlich-Kister solid solution binary model.
    /// The Redlich-Kister model calculates the activity coefficient of the end-members in a
    /// solid solution using the equations:
    /// @f[\ln\gamma_{1}=x_{2}^{2}[a_{0}+a_{1}(3x_{1}-x_{2})+a_{2}(x_{1}-x_{2})(5x_{1}-x_{2})]@f]
    /// and
    /// @f[\ln\gamma_{2}=x_{1}^{2}[a_{0}-a_{1}(3x_{2}-x_{1})+a_{2}(x_{2}-x_{1})(5x_{2}-x_{1})]@f].
    /// The parameters @f$a_0@f$, @f$a_1@f$, and @f$a_2@f$ must be provided.
    /// Set them to zero if not needed.
    /// @param a0 The Redlich-Kister parameter a0
    /// @param a1 The Redlich-Kister parameter a1
    /// @param a2 The Redlich-Kister parameter a2
    auto setChemicalModelRedlichKister(double a0, double a1, double a2) -> MineralPhase&;

    /// Return the MineralMixture instance
    auto mixture() const -> const MineralMixture&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
