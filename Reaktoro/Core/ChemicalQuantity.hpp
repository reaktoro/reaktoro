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
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class ReactionSystem;

/// A class that provides a convenient way to retrieve chemical quantities.
/// Here the term chemical quantity is used in a broad sense. It means any
/// quantity for an element, species, or phase in a chemical system that
/// can be calculated at a chemical state whose temperature, pressure, and
/// mole amounts of all species are known.
///
/// In the example below, the volume of a phase named Gaseous and the pH
/// of the aqueous phase (assuming both phases were defined in the chemical
/// system) are retrieved:
///
/// ~~~
/// ChemicalQuantity q(state);
///
/// const double vol = q["phaseVolume(Gaseous)"];
/// const double pH = q["pH"];
/// ~~~
///
/// The table below shows all possible quantities that can be retrieved from
/// a ChemicalQuantity instance. The first column, **Quantity**, lists the
/// names of the quantities; the second column, **Units**, lists the default
/// units of the quantity; and the third column, **Example**, lists the
/// formatted strings needed to retrieve a quantity.
///
/// | Quantity                 | Units  | Example                                    |
/// | --------                 | -----  | -------                                    |
/// | temperature              | kelvin | `"temperature(units=celsius)"`             |
/// | pressure                 | pascal | `"pressure(units=bar)"`                    |
/// | volume                   | m3     | `"volume(units=cm3)"`                      |
/// | activity                 | ---    | `"activity(CO2(aq))"`                      |
/// | activityCoefficient      | ---    | `"activityCoefficient(Na+)"`               |
/// | fugacity                 | bar    | `"fugacity(CO2(g))"`                       |
/// | chemicalPotential        | J/mol  | `"chemicalPotential(Cl-)"`                 |
/// | elementAmount            | mol    | `"elementAmount(Ca)"`                      |
/// | elementAmountInPhase     | mol    | `"elementAmountInPhase(Mg Aqueous)"`       |
/// | elementMass              | kg     | `"elementMass(Fe units=g)"`                |
/// | elementMassInPhase       | kg     | `"elementMassInPhase(C Gaseous)"`          |
/// | elementMolality          | molal  | `"elementMolality(Cl units=mmolal)"`       |
/// | elementMolarity          | molar  | `"elementMolarity(K units=mmolar)"`        |
/// | speciesAmount            | mol    | `"speciesAmount(H2O(l))"`                  |
/// | speciesMass              | kg     | `"speciesMass(Calcite units=g)"`           |
/// | speciesMoleFraction      | ---    | `"speciesMoleFraction(HCO3-)"`             |
/// | speciesMolality          | molal  | `"speciesMolality(Ca++)"`                  |
/// | speciesMolarity          | molar  | `"speciesMolarity(Mg++)"`                  |
/// | phaseAmount              | mol    | `"phaseAmount(Aqueous)"`                   |
/// | phaseMass                | kg     | `"phaseMass(Dolomite)"`                    |
/// | phaseVolume              | m3     | `"phaseVolume(Gaseous)"`                   |
/// | pH                       | ---    | `"pH"`                                     |
/// | pE                       | ---    | `"pE"`                                     |
/// | Eh                       | volt   | `"Eh"`                                     |
/// | ionicStrength            | molal  | `"ionicStrength"`                          |
/// | fluidVolume              | m3     | `"fluidVolume(units=liter)"`               |
/// | fluidVolumeFraction      | ---    | `"fluidVolumeFraction"`                    |
/// | solidVolume              | m3     | `"solidVolume(units=mm3)"`                 |
/// | solidVolumeFraction      | ---    | `"solidVolumeFraction"`                    |
/// | reactionRate             | mol/s  | `"reactionRate(Dolomite units=mmol/hour)"` |
/// | reactionEquilibriumIndex | ---    | `"reactionEquilibriumIndex(Quartz)"`       |
/// | tag                      | ---    | `"tag"`                                    |
/// | t                        | s      | `"t(units=minute)"`                        |
/// | time                     | s      | `"time(units=year)"`                       |
/// | progress                 | ---    | `"progress"`                               |
///
class ChemicalQuantity
{
public:
    /// A type to describe a chemical quantity function.
    using Function = std::function<double()>;

    /// Disable the default ChemicalQuantity constructor.
    /// This is to enforce the initialization of ChemicalQuantity
    /// instance with a ChemicalSystem instance.
    ChemicalQuantity() = delete;

    /// Construct a ChemicalQuantity instance from a ChemicalSystem object.
    explicit ChemicalQuantity(const ChemicalSystem& system);

    /// Construct a ChemicalQuantity instance from a ReactionSystem object.
    explicit ChemicalQuantity(const ReactionSystem& reactions);

    /// Construct a ChemicalQuantity instance from a ChemicalState object.
    explicit ChemicalQuantity(const ChemicalState& state);

    /// Destroy this ChemicalQuantity instance.
    virtual ~ChemicalQuantity() = default;

    /// Return the chemical system of the ChemicalQuantity instance.
    auto system() const -> const ChemicalSystem&;

    /// Return the chemical reactions of the ChemicalQuantity instance.
    auto reactions() const -> const ReactionSystem&;

    /// Return the chemical state of the ChemicalQuantity instance.
    auto state() const -> const ChemicalState&;

    /// Return the chemical properties of the ChemicalQuantity instance.
    auto properties() const -> const ChemicalProperties&;

    /// Return the reaction rates of the ChemicalQuantity instance.
    auto rates() const -> const ChemicalVector&;

    /// Return the tag variable of the ChemicalQuantity instance.
    auto tag() const -> double;

    /// Update the state of this ChemicalQuantity instance.
    auto update(const ChemicalState& state) -> ChemicalQuantity&;

    /// Update the state of this ChemicalQuantity instance at the time t.
    auto update(const ChemicalState& state, double t) -> ChemicalQuantity&;

    /// Update the state of this ChemicalQuantity instance using provided properties and at the tag.
    auto update(const ChemicalState& state, const ChemicalProperties& properties, double t) -> ChemicalQuantity&;

    /// Return the value of the quantity given as a formatted string.
    auto value(std::string str) const -> double;

    /// Return a created function that calculates the chemical quantity from a formatted string.
    auto function(std::string str) const -> Function;

    /// Return the value of the quantity given as a formatted string.
    auto operator()(std::string str) const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
