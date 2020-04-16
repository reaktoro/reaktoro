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
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>

namespace Reaktoro {
namespace global {

/// A type used to describe all options related to interpolation.
struct InterpolationOptions
{
    /// The units of the temperature values (default: celsius).
    std::string temperature_units = "celsius";

    /// The units of the pressure values (default: bar).
    std::string pressure_units = "bar";

    /// The temperature points used for interpolation.
    std::vector<double> temperatures;

    /// The pressure points used for interpolation.
    std::vector<double> pressures;
};

/// A type used to describe all options related to exception handling.
struct ExceptionOptions
{
    /// The boolean flag that indicates if temperature bounds should be enforced.
    /// When `true`, temperature values that are out of bounds result in runtime error.
    /// When `false`, either the minimum or maximum allowed temperature value is used.
    bool enforce_temperature_bounds = false;

    /// The boolean flag that indicates if pressure bounds should be enforced.
    /// When `true`, pressure values that are out of bounds result in runtime error.
    /// When `false`, either the minimum or maximum allowed pressure value is used.
    bool enforce_pressure_bounds = false;
};

/// A type used to describe all options related to database management.
struct DatabaseOptions
{
    /// When `true`, causes all species with missing data to be ignored during initialization.
    bool exclude_species_with_missing_data = true;
};

/// A type used to describe all options related to aqueous models.
struct AqueousOptions
{
    /// A type used to describe all options related to HKF model.
    struct HKF
    {
        /// The limit on the ionic strength of the model.
        /// This limit value is used instead of larger ionic strengths are specified.
        double ionic_strength_limit = 6.0;
    };

    /// A type used to describe all options related to Debye-Huckel model.
    struct DebyeHuckel
    {
        /// The limit on the ionic strength of the model.
        /// This limit value is used instead of larger ionic strengths are specified.
        double ionic_strength_limit = 6.0;
    };

    /// The options for the HKF model.
    HKF hkf;

    /// The options for the Debye-Huckel model.
    DebyeHuckel debye_huckel;
};

/// A type used to describe all options related to chemical calculations.
struct GlobalOptions
{
    /// The options used for the equilibrium calculations.
    EquilibriumOptions equilibrium;

    /// The options used for the equilibrium path calculations.
    // EquilibriumPathOptions equilibriumpath;

    /// The options used for the kinetic calculations.
    KineticOptions kinetics;

    /// The options used for the property interpolation.
    InterpolationOptions interpolation;

    /// The options used for the exception handling.
    ExceptionOptions exception;

    /// The options used for managing and inializing databases.
    DatabaseOptions database;

    /// The options used for aqueous properties calculations.
    AqueousOptions aqueous;
};

extern GlobalOptions options;

} // namespace global
} // namespace Reaktoro
