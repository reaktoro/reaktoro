// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <memory>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class ReactionSystem;

/// A class that provides a convenient way to calculate chemical quantities.
/// 
/// ~~~
/// ChemicalQuantity quantity(state);
/// double vol = quantity["phaseVolume(Aqueous)"];
/// ~~~
/// 
/// | Quantity | Units | Example | 
/// | -------- | ----- | ------- |
/// | eh | volt | `"eH"` |
/// | elementamount | mol | `"elementAmount(Na units=mmol)"` |
/// | elementamountinphase | mol | `"elementAmountInPhase(C Aqueous)"` |
/// | elementmass | kg | `"elementMass(Ca units=g)"` |
/// | elementmassinphase | kg | `"elementMassInPhase(C Gaseous)"` |
/// | elementmolality | molal | `"elementMolality(C)"` |
/// | elementmolarity | molar | `"elementMolarity(C)"` |
/// | fluidvolume | m3 | `"fluidVolume"` |
/// | fugacity | bar | `"fugacity(CO2(g))"` |
/// | ionicstrength | molal | `"ionicStrength(units=mmolal)"` |
/// | phaseamount | mol | `"phaseAmount(Calcite)"` |
/// | phasemass | kg | `"phaseMass(Quartz units=mg)"` |
/// | phasevolume | m3 | `"phaseVolume(Aqueous units=liter)"` |
/// | pressure | pascal | `"pressure(units=bar)"` |
/// | reactionrate | mol/s | `"reactionRate(units=mmol/day)"` |
/// | solidvolume | m3 | `"solidVolume(units=cm3)"` |
/// | speciesamount | mol | `"speciesAmount(H2O(l))"` |
/// | speciesmass | kg | `"speciesMass(CO2(g) units=mg)"` |
/// | speciesmolality | molal | `"speciesMolality(HCO3- units=mmolal)"` |
/// | speciesmolarity | molar | `"speciesMolarity(Cl-)"` |
/// | t | s | `"t(units=year)"` |
/// | temperature | kelvin | `"temperature(units=celsius)"` |
/// | time | s | `"time(units=minute)"` |
/// | volume | m3 | `"volume(units=mm3)"` |
class ChemicalQuantity
{
public:
    /// A type to describe a chemical quantity function.
    using Function = std::function<double()>;

    /// Construct a default ChemicalQuantity instance.
    ChemicalQuantity();

    /// Construct a copy of a ChemicalQuantity instance.
    ChemicalQuantity(const ChemicalQuantity& other);

    /// Construct a ChemicalQuantity instance from a ChemicalSystem object.
    explicit ChemicalQuantity(const ChemicalSystem& system);

    /// Construct a ChemicalQuantity instance from a ReactionSystem object.
    explicit ChemicalQuantity(const ReactionSystem& reactions);

    /// Destroy this ChemicalQuantity instance.
    virtual ~ChemicalQuantity();

    /// Assign a ChemicalQuantity instance to this.
    auto operator=(ChemicalQuantity other) -> ChemicalQuantity&;

    /// Update the state of this ChemicalQuantity instance.
    auto update(const ChemicalState& state) -> void;

    /// Update the state of this ChemicalQuantity instance.
    auto update(const ChemicalState& state, double t) -> void;

    /// Return the value of the quantity given as a formatted string.
    auto value(std::string str) const -> double;

    /// Return a created function that calculates the chemical quantity from a formatted string.
    auto function(std::string str) const -> Function;

    /// Return the value of the quantity given as a formatted string.
    auto operator[](std::string str) const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
