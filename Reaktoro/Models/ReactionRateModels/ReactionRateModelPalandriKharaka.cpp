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

#include "ReactionRateModelPalandriKharaka.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/SurfaceList.hpp>
#include <Reaktoro/Serialization/Models/ReactionRateModels.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {
namespace detail {

using Catalyst = ReactionRateModelParamsPalandriKharaka::Catalyst;
using Mechanism = ReactionRateModelParamsPalandriKharaka::Mechanism;

/// Construct a function that computes the activity-based contribution of a catalyst in the mineral reaction rate.
auto mineralCatalystFnActivity(Catalyst const& catalyst, ReactionRateModelGeneratorArgs args) -> Fn<real(ChemicalProps const&)>
{
    auto const& formula = catalyst.formula;
    auto const& power = catalyst.power;
    auto const& species = args.species;

    auto const aqspecies = species.withAggregateState(AggregateState::Aqueous);
    auto const iaqueousspecies = aqspecies.findWithFormula(formula);

    if(aqspecies.size() == 0 || iaqueousspecies >= aqspecies.size())
    {
        // warningif(true, "Ignoring Palandri-Kharaka catalytic effect (based on activity) in mineral reaction rate because no aqueous species with formula `", formula, "` exists in the aqueous phase of the system.");
        return [](ChemicalProps const& props) { return 1.0; };
    }

    auto const name = aqspecies[iaqueousspecies].name();
    auto const ispecies = species.findWithName(name);

    auto fn = [=](ChemicalProps const& props)
    {
        auto const& ai = props.speciesActivity(ispecies);
        return pow(ai, power);
    };

    return fn;
}

/// Construct a function that computes the partial-pressure-based contribution of a catalyst in the mineral reaction rate.
auto mineralCatalystFnPartialPressure(Catalyst const& catalyst, ReactionRateModelGeneratorArgs args) -> Fn<real(ChemicalProps const&)>
{
    auto const& formula = catalyst.formula;
    auto const& power = catalyst.power;
    auto const& species = args.species;

    auto const gases = species.withAggregateState(AggregateState::Gas);
    auto const igas = gases.findWithFormula(formula);

    if(gases.size() == 0 || igas >= gases.size())
    {
        // warningif(true, "Ignoring Palandri-Kharaka catalytic effect (based on partial pressure) in mineral reaction rate because no gaseous species with formula `", formula, "` exists in the gaseous phase of the system.");
        return [](ChemicalProps const& props) { return 1.0; };
    }

    auto const name = gases[igas].name();
    auto const ispecies = species.findWithName(name);

    auto fn = [=](ChemicalProps const& props)
    {
        auto const P  = props.pressure(); // pressure in Pa
        auto const xi = props.speciesMoleFraction(ispecies);
        auto const Pi = xi * P * 1e-5; // partial pressure in bar!
        return pow(Pi, power);
    };

    return fn;
}

/// Construct a function that computes the contribution of a catalyst in the mineral reaction rate.
auto mineralCatalystFn(Catalyst const& catalyst, ReactionRateModelGeneratorArgs args) -> Fn<real(ChemicalProps const&)>
{
    if(catalyst.property == "a")
        return mineralCatalystFnActivity(catalyst, args);
    if(catalyst.property == "P")
        return mineralCatalystFnPartialPressure(catalyst, args);
    errorif(true, "Expecting mineral catalyst property symbol to be either `a` or `P`, but got `", catalyst.property, "` instead.");
}

auto mineralMechanismFn(Mechanism const& mechanism, ReactionRateModelGeneratorArgs args) -> Fn<real(ChemicalProps const&)>
{
    // The universal gas constant (in J/(mol*K))
    const auto R = universalGasConstant;

    // Create the mineral catalyst functions
    Vec<Fn<real(ChemicalProps const&)>> catalyst_fns;
    for(auto&& catalyst : mechanism.catalysts)
        catalyst_fns.push_back(mineralCatalystFn(catalyst, args));

    // The name of the mineral from the name of the reaction
    const auto mineral = args.name;

    // Define the mineral mechanism function
    auto fn = [=](ChemicalProps const& props)
    {
        const auto& aprops = AqueousProps::compute(props);

        const auto& lgk = mechanism.lgk.value();
        const auto& E = mechanism.E.value() * 1e3; // from kJ to J
        const auto& p = mechanism.p.value();
        const auto& q = mechanism.q.value();

        const auto T = props.temperature();
        const auto k0 = pow(10, lgk);
        const auto k = k0 * exp(-E/R * (1.0/T - 1.0/298.15));

        const auto Omega  = aprops.saturationRatio(mineral); // TODO: Find a way to determine the index of the mineral once in aprops.saturationSpecies to avoid index search in every call of Palandri-Kharaka model (overhead here should be minimal though compared to other operations).
        const auto pOmega = p != 1.0 ? pow(Omega, p) : Omega;
        const auto qOmega = q != 1.0 ? pow(1 - pOmega, q) : 1 - pOmega;

        real g = 1.0;
        for(auto const& catalystfn : catalyst_fns)
            g *= catalystfn(props);

        return k * qOmega * g;
    };

    return fn;
}

} // namespace detail

auto ReactionRateModelPalandriKharaka() -> ReactionRateModelGenerator
{
    const auto params = Params::embedded("PalandriKharaka.json");
    return ReactionRateModelPalandriKharaka(params);
}

auto ReactionRateModelPalandriKharaka(Params const& params) -> ReactionRateModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ReactionRateModelParams"), "Expecting Palandri-Kharaka mineral rate parameters in given Params object, but it lacks a `ReactionRateModelParams` section within which another section `PalandriKharaka` should exist.");
    errorif(!data.at("ReactionRateModelParams").exists("PalandriKharaka"), "Expecting Palandri-Kharaka mineral rate parameters in given Params object, under the section `PalandriKharaka`.");
    errorif(!data.at("ReactionRateModelParams").at("PalandriKharaka").isDict(), "Expecting section `PalandriKharaka` with Palandri-Kharaka mineral rate parameters to be a dictionary.");

    Vec<ReactionRateModelParamsPalandriKharaka> paramsvec;
    for(auto const& [key, value] : data["ReactionRateModelParams"]["PalandriKharaka"].asDict())
        paramsvec.push_back(value.as<ReactionRateModelParamsPalandriKharaka>());

    return ReactionRateModelPalandriKharaka(paramsvec);
}

auto ReactionRateModelPalandriKharaka(ReactionRateModelParamsPalandriKharaka const& params) -> ReactionRateModelGenerator
{
    ReactionRateModelGenerator model = [=](ReactionRateModelGeneratorArgs args)
    {
        Vec<Fn<real(ChemicalProps const&)>> mechanism_fns;
        for(auto const& mechanism : params.mechanisms)
            mechanism_fns.push_back(detail::mineralMechanismFn(mechanism, args));

        const auto imineralsurface = args.surfaces.indexWithName(args.name);

        ReactionRateModel fn = [=](ChemicalProps const& props) -> ReactionRate
        {
            const auto area = props.surfaceArea(imineralsurface);
            real sum = 0.0;
            for(auto const& mechanismfn : mechanism_fns)
                sum += mechanismfn(props);
            return area * sum;
        };

        return fn;
    };

    return model;
}

auto ReactionRateModelPalandriKharaka(Vec<ReactionRateModelParamsPalandriKharaka> const& paramsvec) -> ReactionRateModelGenerator
{
    ReactionRateModelGenerator model = [=](ReactionRateModelGeneratorArgs args) -> ReactionRateModel
    {
        const auto mineral = args.name;
        const auto idx = indexfn(paramsvec, RKT_LAMBDA(x, x.mineral == mineral || contains(x.othernames, mineral)));
        errorif(idx >= paramsvec.size(), "Could not find a mineral with name `", mineral, "` in the provided set of Palandri-Kharaka parameters.");
        const auto params = paramsvec[idx];
        return ReactionRateModelPalandriKharaka(params)(args);
    };

    return model;
}

} // namespace Reaktoro
