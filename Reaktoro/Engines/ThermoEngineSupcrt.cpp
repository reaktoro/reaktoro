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

#include "ThermoEngineSupcrt.hpp"

// C++ includes
#include <any>

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/OptimizationUtils.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Databases/DatabaseSupcrt.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesElectroState.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>

namespace Reaktoro {
namespace {

/// The signature of a function that calculates the thermodynamic properties of water
using WaterThermoPropsFn = std::function<WaterThermoState(double, double)>;

/// The signature of a function that calculates the electrostatic properties of water
using WaterElectroPropsFn = std::function<WaterElectroState(double, double)>;

/// The signature of a function that calculates the standard thermodynamic properties of a species
using SpeciesThermoPropsFn = std::function<SpeciesThermoState(double, double, std::string)>;

/// Convert a SpeciesThermoState object into a StandardThermoProps one
auto convert(const SpeciesThermoState& state) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = state.gibbs_energy;
    props.A0  = state.helmholtz_energy;
    props.U0  = state.internal_energy;
    props.H0  = state.enthalpy;
    props.S0  = state.entropy;
    props.V0  = state.volume;
    props.Cp0 = state.heat_capacity_cp;
    props.Cv0 = state.heat_capacity_cv;
    return props;
}

/// The standard thermodynamic model function object based on SUPCRT
struct SupcrtStandardThermoModelFn
{
    /// The Haar--Gallagher--Kell (1984) equation of state for water
    WaterThermoPropsFn water_thermo_props_hgk_fn;

    /// The Wagner and Pruss (1995) equation of state for water
    WaterThermoPropsFn water_thermo_props_wagner_pruss_fn;

    /// The Johnson and Norton equation of state for the electrostatic state of water
    WaterElectroPropsFn water_eletro_props_fn;

    /// The HKF equation of state for the thermodynamic state of aqueous, gaseous and mineral species
    SpeciesThermoPropsFn species_thermo_props_hkf_fn;

    /// Construct a SupcrtStandardThermoModelFn object
    SupcrtStandardThermoModelFn()
    {
        // Initialize the Haar--Gallagher--Kell (1984) equation of state for water
        water_thermo_props_hgk_fn = [](real T, real P)
        {
            return Reaktoro::waterThermoStateHGK(T, P, StateOfMatter::Liquid);
        };

        water_thermo_props_hgk_fn = memoize(water_thermo_props_hgk_fn);

        // Initialize the Wagner and Pruss (1995) equation of state for water
        water_thermo_props_wagner_pruss_fn = [](real T, real P)
        {
            return Reaktoro::waterThermoStateWagnerPruss(T, P, StateOfMatter::Liquid);
        };

        water_thermo_props_wagner_pruss_fn = memoize(water_thermo_props_wagner_pruss_fn);

        // Initialize the Johnson and Norton equation of state for the electrostatic state of water
        water_eletro_props_fn = [=](double T, double P)
        {
            const WaterThermoState wts = water_thermo_props_wagner_pruss_fn(T, P);
            return waterElectroStateJohnsonNorton(T, P, wts);
        };

        water_eletro_props_fn = memoize(water_eletro_props_fn);
    }

    /// Return the standard thermodynamic properties of a species at given temperature and pressure.
    auto operator()(real T, real P, const Species& species) const -> StandardThermoProps
    {
        const auto& data = species.attachedData();
        if(data.type() == typeid(ParamsMaierKelly))
            return standardThermoPropsMaierKelly(T, P, species);
        if(data.type() == typeid(ParamsMaierKellyHKF))
            return standardThermoPropsMaierKellyHKF(T, P, species);
        if(isAlternativeWaterName(species.name()))
            return standardThermoPropsAqueousSolventHKF(T, P, species);
        if(data.type() == typeid(ParamsAqueousSoluteHKF))
            return standardThermoPropsAqueousSoluteHKF(T, P, species);
        RuntimeError("Failure at SupcrtStandardThermoModelFn::operator().",
            "There is no SUPCRT parameters for species named " + species.name() + ".");
        return {};
    }

    /// Return the standard thermodynamic properties of the aqueous solvent water using HKF model.
    auto standardThermoPropsAqueousSolventHKF(real T, real P, const Species& species) const -> StandardThermoProps
    {
        const auto wts = water_thermo_props_wagner_pruss_fn(T, P);
        const auto res = speciesThermoStateSolventHKF(T, P, wts);
        return convert(res);
    }

    /// Return the standard thermodynamic properties of an aqueous solute using HKF model.
    auto standardThermoPropsAqueousSoluteHKF(real T, real P, const Species& species) const -> StandardThermoProps
    {
        const auto params = std::any_cast<ParamsAqueousSoluteHKF>(species.attachedData());
        const auto wts = water_thermo_props_wagner_pruss_fn(T, P);
        const auto wes = water_eletro_props_fn(T, P);
        const auto g = functionG(T, P, wts);
        const auto aes = speciesElectroStateHKF(g, params);
        const auto res = speciesThermoStateSoluteHKF(T, P, params, aes, wes);
        return convert(res);
    }

    /// Return the standard thermodynamic properties of a species using Maier-Kelly model.
    auto standardThermoPropsMaierKelly(real T, real P, const Species& species) const -> StandardThermoProps
    {
        const auto params = std::any_cast<ParamsMaierKelly>(species.attachedData());
        const auto res = speciesThermoStateHKF(T, P, params);
        return convert(res);
    }

    /// Return the standard thermodynamic properties of a mineral species using Maier-Kelly-HKF model.
    auto standardThermoPropsMaierKellyHKF(real T, real P, const Species& species) const -> StandardThermoProps
    {
        const auto params = std::any_cast<ParamsMaierKellyHKF>(species.attachedData());
        const auto res = speciesThermoStateHKF(T, P, params);
        return convert(res);
    }
};

} // namespace

ThermoEngineSupcrt::ThermoEngineSupcrt(const DatabaseSupcrt& database)
: ThermoEngine(database, SupcrtStandardThermoModelFn())
{
}

} // namespace Reaktoro
