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
#include <map>
#include <memory>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species from SUPCRT databases.
/// @see Element, Species
/// @ingroup Databases
class DatabaseSupcrt : public Database
{
public:
    /// Construct a default DatabaseSupcrt object.
    DatabaseSupcrt();

    /// Construct a DatabaseSupcrt object with given name of a built-in database file.
    /// If `name` does not correspond to one of the following names, an exception is thrown:
    /// - `supcrt98`
    /// - `supcrt07`
    /// - `supcrt98-organics`
    /// - `supcrt07-organics`
    /// @param name The name of the built-in SUPCRT database.
    static auto withName(std::string name) ->  DatabaseSupcrt;

    /// Construct a DatabaseSupcrt object with given database file.
    /// If `path` does not point to a valid database file, an exception is thrown.
    /// @param path The path, including file name, to the database file.
    static auto fromFile(std::string path) ->  DatabaseSupcrt;
};

/// The type used to store HKF parameters based on SUPCRT databases for an aqueous solute.
struct SupcrtParamsAqueousSoluteHKF
{
    /// The electrical charge of the aqueous solute
    double charge = 0.0;

    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf = 0.0;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf = 0.0;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    double Sr = 0.0;

    /// The coefficient a1 of the HKF equation of state of the aqueous solute (in units of cal/(mol*bar))
    double a1 = 0.0;

    /// The coefficient a2 of the HKF equation of state of the aqueous solute (in units of cal/mol)
    double a2 = 0.0;

    /// The coefficient a3 of the HKF equation of state of the aqueous solute (in units of (cal*K)/(mol*bar))
    double a3 = 0.0;

    /// The coefficient a4 of the HKF equation of state of the aqueous solute (in units of (cal*K)/mol)
    double a4 = 0.0;

    /// The coefficient c1 of the HKF equation of state of the aqueous solute (in units of cal/(mol*K))
    double c1 = 0.0;

    /// The coefficient c2 of the HKF equation of state of the aqueous solute (in units of (cal*K)/mol)
    double c2 = 0.0;

    /// The conventional Born coefficient of the aqueous solute at reference temperature 298.15 K and pressure 1 bar (in units of cal/mol)
    double wref = 0.0;

    /// The dissociation of a neutral aqueous solute into charged species.
    std::map<std::string, double> dissociation; // TODO: This dissociation data should not exist here, but in a more central location that will work for other databases as well.
};

/// The type used to store Maier-Kelly parameters based on SUPCRT databases for gaseous/liquid species.
struct SupcrtParamsMaierKelly
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf = 0.0;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf = 0.0;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    double Sr = 0.0;

    /// The coefficient a of the Maier-Kelly equation of state of the gaseous/liquid species (in units of cal/(mol*K))
    double a = 0.0;

    /// The coefficient b of the Maier-Kelly equation of state of the gaseous/liquid species (in units of cal/(mol*K^2))
    double b = 0.0;

    /// The coefficient c of the Maier-Kelly equation of state of the gaseous/liquid species (in units of (cal*K)/mol)
    double c = 0.0;

    /// The maximum temperature at which the Maier-Kelly equation of state can be applied for the gaseous species (in units of K)
    double Tmax = 0.0;

    // The critical temperature of the species (in units of K)
    double Tcr = 0.0; // TODO: All these critical properties should not exist here, but in a more central location that will work for other databases as well.

    // The critical pressure of the species (in units of Pa)
    double Pcr = 0.0;

    // The acentric factor of the species
    double acentric_factor = 0.0;
};

/// The type used to store Maier-Kelly-HKF parameters based on SUPCRT databases for mineral species.
struct SupcrtParamsMaierKellyHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf = 0.0;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf = 0.0;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    double Sr = 0.0;

    /// The standard molal volume of the mineral species at reference temperature and pressure (in units of cm^3/mol)
    double Vr = 0.0;

    /// The number of phase transitions that the mineral may undergo
    int nptrans = 0;

    /// The coefficients ai of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of cal/(mol*K))
    std::vector<double> a;

    /// The coefficients bi of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of cal/(mol*K^2))
    std::vector<double> b;

    /// The coefficients ci of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of (cal*K)/mol)
    std::vector<double> c;

    /// The temperatures at which the mineral experiences phase transition along the line of reference pressure (in units of K)
    std::vector<double> Ttr;

    /// The change in the standard molal enthalpy of each mineral phase transition (in units of cal/mol)
    std::vector<double> Htr;

    /// The change in the standard molal volume of each mineral phase transition (in units of cm^3/mol)
    std::vector<double> Vtr;

    /// The Clapeyron slote at each mineral phase transition (in units of bar/K)
    std::vector<double> dPdTtr;

    /// The maximum temperature at which the Maier-Kelly-HKF equation of state can be applied for the mineral (in units of K)
    double Tmax = 0.0;
};

} // namespace Reaktoro
