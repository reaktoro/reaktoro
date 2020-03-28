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

#include "ThermoEngineThermoFun.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Databases/DatabaseThermoFun.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Convert ThermoFun::real to Reaktoro::real
auto convertScalar(Reaktoro_::real other) -> real
{
    real converted;
    converted.val = other.val;
    converted.ddP = other.ddp;
    converted.ddT = other.ddt;
    return converted;
}

/// Convert a SpeciesThermoState object into a StandardThermoProps one
auto convertProps(const ThermoFun::ThermoPropertiesSubstance& other) -> StandardThermoProps
{
    StandardThermoProps converted;
    converted.G0  = convertScalar(other.gibbs_energy);
    converted.A0  = convertScalar(other.helmholtz_energy);
    converted.U0  = convertScalar(other.internal_energy);
    converted.H0  = convertScalar(other.enthalpy);
    converted.S0  = convertScalar(other.entropy);
    converted.V0  = convertScalar(other.volume*1e-05); // from J/bar to m3/mol
    converted.Cp0 = convertScalar(other.heat_capacity_cp);
    converted.Cv0 = convertScalar(other.heat_capacity_cv);
    return converted;
}

/// The standard thermodynamic model function object based on ThermoFun
struct ThermoFunStandardThermoModelFn
{
    // The ThermoFun::ThermoEngine object
    mutable ThermoFun::ThermoEngine thermofun_engine; // mutable needed because some methods in ThermoFun is not const

    /// Construct a ThermoFunStandardThermoModelFn object
    ThermoFunStandardThermoModelFn(const DatabaseThermoFun& database)
    : thermofun_engine(std::any_cast<ThermoFun::Database>(database.data()))
    {
        // Set solvent symbol, the HGK, JN water solvent model are defined in this record
        thermofun_engine.setSolventSymbol("H2O@");
    }

    /// Return the standard thermodynamic properties of a species at given temperature and pressure.
    auto operator()(Temperature T, Pressure P, const Species& species) const -> StandardThermoProps
    {
        double Tval = T;
        double Pval = P;
        const auto& db = thermofun_engine.database();
        if(db.containsSubstance(species.name()))
        {
            const auto props = thermofun_engine.thermoPropertiesSubstance(Tval, Pval, species.name());
            return convertProps(props);
        }
        RuntimeError("Failure at ThermoFunStandardThermoModelFn::operator().",
            "There is no species named " + species.name() + " in the ThermoFun database.");
        return {};
    }
};

} // namespace

ThermoEngineThermoFun::ThermoEngineThermoFun(const DatabaseThermoFun& database)
: ThermoEngine(database, ThermoFunStandardThermoModelFn(database))
{
}

} // namespace Reaktoro
