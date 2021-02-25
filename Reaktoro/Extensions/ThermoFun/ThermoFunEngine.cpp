// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "ThermoFunEngine.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Convert a SpeciesThermoState object into a StandardThermoProps one
auto convertProps(const ThermoFun::ThermoPropertiesSubstance& other) -> StandardThermoProps
{
    StandardThermoProps converted;
    converted.G0  = other.gibbs_energy.val;
    converted.H0  = other.enthalpy.val;
    converted.V0  = other.volume.val * 1.0e-05; // from J/bar to m3/mol
    converted.Cp0 = other.heat_capacity_cp.val;
    converted.Cv0 = other.heat_capacity_cv.val;
    return converted;
}

} // namespace

struct ThermoFunEngine::Impl
{
    /// The ThermoFun::ThermoEngine object.
    mutable ThermoFun::ThermoEngine engine; // mutable needed because some methods in ThermoFun is not const

    /// The ThermoFun::Database object.
    ThermoFun::Database database;

    /// Costruct a Impl object with given ThermoFun::Database object.
    Impl(const ThermoFun::Database& database)
    : engine(database), database(database)
    {
        // Set solvent symbol, the HGK, JN water solvent model are defined in this record
        engine.setSolventSymbol("H2O@");
    }

    /// Return the standard thermodynamic properties of a chemical species with given name.
    auto props(double T, double P, const String& species) const -> StandardThermoProps
    {
        if(database.containsSubstance(species))
        {
            const auto props = engine.thermoPropertiesSubstance(T, P, species);
            return convertProps(props);
        }
        error("Cannot evaluate the standard thermodynamic properties of species with "
            "name ", species, " because it is know available in the ThermoFun database.");
        return {};
    }
};


ThermoFunEngine::ThermoFunEngine(const ThermoFun::Database& database)
: pimpl(new Impl(database))
{}

auto ThermoFunEngine::database() const -> const ThermoFun::Database&
{
    return pimpl->database;
}

auto ThermoFunEngine::props(real T, real P, const String& species) const -> StandardThermoProps
{
    return pimpl->props(T, P, species);
}

} // namespace Reaktoro
