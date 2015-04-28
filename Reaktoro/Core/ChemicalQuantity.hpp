// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class ReactionSystem;

/// A class that facilitates the calculation of chemical quantities using formatted strings.
/// The table below shows the chemical quantities that can be calculated with
/// its respective formatted string:
///
/// | String           | Quantity                                | Example                  |
/// |:----------------:|:---------------------------------------:|:------------------------:|
/// | `n[`*species*`]` | the molar amount of species named *species* | `n[H2O(l)]`, `n[CO2(g)]` |
/// | `b[`*element*`]` | the molar amount of element named *element* | `b[H]`, `b[C]`, `b[O]`   |
/// | `b[`*element*`][`*phase*`]` | the molar amount of element named *element* in a phase named *phase* | `b[Na][Aqueous]`, `b[C][Gaseous]`|
/// | `x[`*species*`]` | the molar fraction of species named *species* | `x[H2O(l)]`, `x[CO2(g)]` |
/// | `m[`*species*`]` | the molality of the aqueous species named *species* | `m[Na+]`, `m[HCO3-]` |
/// | `m[`*element*`]` | the molality of the element *element* in the aqueous phase | `m[Na]`, `m[C]` |
/// | `a[`*species*`]` | the activity of the species named *species* | `a[H+]`, `a[CO2(aq)]` |
/// | `g[`*species*`]` | the activity coefficient of the species named *species* | `g[H+]`, `g[CO2(aq)]` |
/// | `pH` | the pH of the aqueous phase | `pH` |
///
/// The above strings can be combined with the units of the extracted quantity. For example, the molar amount of species `H+` can be extracted in `mmol` units as `n[H+]:mmol`.
///
/// Note that extracting the pH of the aqeous phase and the molality of one of its species will only succeed as long as the aqueous phase is called `Aqueous`.
///
/// **Usage:**
/// ~~~
/// ChemicalQuantity quantity(system);
/// quantity.update(state);
/// std::cout << "n[CO2(aq)] = " << quantity.value("n[CO2(aq)]") << std::endl;
/// std::cout << "pH = " << quantity["pH"] << std::endl;
/// ~~~
class ChemicalQuantity
{
public:
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
    auto value(std::string quantity) const -> double;

    /// Return the value of the quantity given as a formatted string.
    auto operator[](std::string quantity) const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
