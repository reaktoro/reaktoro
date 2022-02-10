// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

auto convertProps(const ThermoFun::ThermoPropertiesSubstance& other, const ThermoFun::Substance& substance) -> StandardThermoProps
{
    StandardThermoProps converted;
    converted.G0  = other.gibbs_energy.val;
    converted.H0  = other.enthalpy.val;
    converted.V0  = other.volume.val * 1.0e-05; // from J/bar to m3/mol
    converted.Cp0 = other.heat_capacity_cp.val;

    // Reaktoro expects zero standard molar volumes for gases.
    if(substance.aggregateState() == ThermoFun::AggregateState::GAS)
        converted.V0 = 0.0;

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
    auto props(const real& T, const real& P, const String& species) const -> StandardThermoProps
    {
        auto const& substances = database.mapSubstances();
        const auto it = substances.find(species);

        errorif(it == substances.end(), "Expecting a species name that exists in the ThermoFun database, but got `", species, "` instead.");

        const auto substance = it->second;

        return props(T, P, substance);
    }

    /// Return the standard thermodynamic properties of a chemical species with given ThermoFun::Substance object.
    auto props(const real& T, const real& P, const ThermoFun::Substance& substance) const -> StandardThermoProps
    {
        double Tval = T.val();
        double Pval = P.val();
        const auto props = engine.thermoPropertiesSubstance(Tval, Pval, substance);
        return convertProps(props, substance);
    }
};

ThermoFunEngine::ThermoFunEngine(const ThermoFun::Database& database)
: pimpl(new Impl(database))
{}

auto ThermoFunEngine::database() const -> const ThermoFun::Database&
{
    return pimpl->database;
}

auto ThermoFunEngine::props(const real& T, const real& P, const String& species) const -> StandardThermoProps
{
    return pimpl->props(T, P, species);
}

auto ThermoFunEngine::props(const real& T, const real& P, const ThermoFun::Substance& substance) const -> StandardThermoProps
{
    return pimpl->props(T, P, substance);
}

} // namespace Reaktoro
