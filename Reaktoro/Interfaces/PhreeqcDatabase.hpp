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
#include <set>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>

namespace Reaktoro {

// Forward declarations
class Database;
class Element;
class Species;

class PhreeqcDatabase
{
public:
    /// Construct a default PhreeqcDatabase instance
    PhreeqcDatabase();

    /// Construct a custom PhreeqcDatabase instance
    /// @param filename The path to the Phreeqc database file
    explicit PhreeqcDatabase(std::string filename);

    /// Load a Phreeqc database.
    /// @param filename The path to the Phreeqc database file
    auto load(std::string filename) -> void;

    auto numElements() const -> unsigned;

    auto numAqueousSpecies() const -> unsigned;

    auto numGaseousSpecies() const -> unsigned;

    auto numMineralSpecies() const -> unsigned;

    auto numMasterSpecies() const -> unsigned;

    auto numProductSpecies() const -> unsigned;

    auto element(Index index) const -> Element;

    auto elements() const -> const std::vector<Element>&;

    auto aqueousSpecies(Index index) const -> Species;

    auto aqueousSpecies(std::string name) const -> Species;

    auto aqueousSpecies() const -> const std::vector<Species>&;

    auto gaseousSpecies(Index index) const -> Species;

    auto gaseousSpecies(std::string name) const -> Species;

    auto gaseousSpecies() const -> const std::vector<Species>&;

    auto mineralSpecies(Index index) const -> Species;

    auto mineralSpecies(std::string name) const -> Species;

    auto containsAqueousSpecies(std::string name) const -> bool;

    auto containsGaseousSpecies(std::string name) const -> bool;

    auto containsMineralSpecies(std::string name) const -> bool;

    auto mineralSpecies() const -> const std::vector<Species>&;

    auto masterSpecies() const -> std::set<std::string>;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// A type for storing Phreeqc parameters of a species.
struct SpeciesThermoParamsPhreeqc
{
    struct ReactionParams
    {
        /// The reaction equation defining the Phreeqc species in terms of master species.
        ReactionEquation equation;

        /// The equilibrium constant of the product species at 25 °C.
        double log_k;

        /// The standard enthalpy of the reaction at 25 °C (in units of kJ/mol).
        /// This parameter is used in the Van't Hoff equation to calculate the equilibrium constant
        /// of the reaction at temperature @f$T@f$:
        /// @f[\ln K=\ln K^{298.15\mathrm{K}}-\frac{\Delta H^{\circ}}{R}\left(\frac{1}{T}-\frac{1}{298.15}\right)@f],
        /// where @f$R@f$ is the universal gas constant. This equation requires the standard enthalpy of reaction
        /// and its equilibrium constant at 25 °C.
        double delta_h;

        /// The coefficients of the analytical expression of the equilibrium constant of the reaction.
        /// The analytical expression is:
        /// @f[\log_{10}K=A_{1}+A_{2}T+\frac{A_{3}}{T}+A_{4}\log_{10}T+\frac{A_{5}}{T^{2}}+A_{6}T^{2}@f],
        /// where @f$T@f$ is temperature in Kelvin, and @f$A_i@f$ are the coefficients.
        std::vector<double> analytic;
    };

    ReactionParams reaction;
};

} // namespace Reaktoro
