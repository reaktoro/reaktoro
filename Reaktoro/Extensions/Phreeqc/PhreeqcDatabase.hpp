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
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class that reads and retrieves data from a PHREEQC thermodynamic database.
class PhreeqcDatabase
{
public:
    /// Construct a default PhreeqcDatabase object.
    PhreeqcDatabase();

    /// Construct a PhreeqcDatabase object with given database. This
    /// constructor supports initialization of the PhreeqcDatabase object with
    /// either a path to the database file, including its file name, or a
    /// multi-line string containing the database contents itself.
    /// @param database The path to the database file or its contents in a string
    explicit PhreeqcDatabase(String database);

    /// Return the elements in the database.
    auto elements() const -> ElementListConstRef;

    /// Return the species in the database.
    auto species() const -> SpeciesListConstRef;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// The data of a chemical species in any PHREEQC database.
struct PhreeqcSpeciesData
{
    /// The reactants and their stoichiometric coefficientsformation reaction of the species.
    Pairs<String, double> reactants;

    /// The equilibrium constant of the formation reaction at 25 °C.
    real log_k = 0.0;

    /// The standard enthalpy of the formation reaction at 25 °C (in kJ/mol).
    /// This parameter is used in the Van't Hoff equation to calculate the
    /// equilibrium constant of the reaction at temperature *T*:
    /// @f[\ln K=\ln K^{298.15\mathrm{K}}-\frac{\Delta H^{\circ}}{R}\left(\frac{1}{T}-\frac{1}{298.15}\right)@f],
    /// where *R* is the universal gas constant. This equation requires the
    /// standard enthalpy of reaction and its equilibrium constant at 25 °C.
    double delta_h;

    /// The coefficients of equilibrium constant function of the formation reaction.
    /// The analytical expression is:
    /// @eqc{\log_{10}K=A_{1}+A_{2}T+\frac{A_{3}}{T}+A_{4}\log_{10}T+\frac{A_{5}}{T^{2}}+A_{6}T^{2}},
    /// where @f$T@f$ is temperature in K, and @f$A_i@f$ are the coefficients.
    Vec<double> analytic;
};

} // namespace Reaktoro
