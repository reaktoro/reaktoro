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
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ScalarTypes.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class Partition;
class ReactionSystem;

/// A type that defines a kinetic problem
class KineticProblem
{
public:
    /// Construct a default KineticProblem instance
    KineticProblem();

    /// Construct a KineticProblem instance with given reactions.
    /// By default, all species that do not participate in the system of reactions
    /// are assumed to be inert. The method KineticProblem::setPartition offers more
    /// control on the groupt of kinetic, equilibrium, and inert species.
    /// @param reactions The reactions represented as a ReactionSystem instance
    explicit KineticProblem(const ReactionSystem& reactions);

    /// Construct a copy of a KineticProblem instance
    KineticProblem(const KineticProblem& other);

    /// Destroy this KineticProblem instance
    virtual ~KineticProblem();

    /// Assign a KineticProblem instance to this
    auto operator=(KineticProblem other) -> KineticProblem&;

    /// Set the temperature for the kinetic calculation (in units of K).
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value (in units of K)
    auto setTemperature(double val) -> KineticProblem&;

    /// Set the temperature for the kinetic calculation with given units.
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value
    /// @param units The units of the temperature (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine)
    auto setTemperature(double val, std::string units) -> KineticProblem&;

    /// Set the pressure for the kinetic calculation (in units of Pa)
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value (in units of Pa)
    /// @param units The units of the pressure (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine)
    auto setPressure(double val) -> KineticProblem&;

    /// Set the pressure for the kinetic calculation (in units of Pa)
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value
    /// @param units The units of the pressure (Pa, kPa, MPa, GPa, atm, mmHg, inHg, psi, kpsi, Mpsi, psf, bar, torr, inH2O, ftH2O, pascal)
    auto setPressure(double val, std::string units) -> KineticProblem&;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(const Partition& partition) -> KineticProblem&;

    /// Return the temperature for the kinetic calculation (in units of K)
    auto temperature() const -> double;

    /// Return the pressure for the kinetic calculation (in units of Pa)
    auto pressure() const -> double;

    /// Return a reference to the ReactionSystem instance used to create this KineticProblem instance
    auto reactions() const -> const ReactionSystem&;

    /// Return a reference to the ChemicalSystem instance used to create this KineticProblem instance
    auto system() const -> const ChemicalSystem&;

    /// Return a reference to the Partition instance used to create this KineticProblem instance
    auto partition() const -> const Partition&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
