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

#include "AqueousPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDrummondCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDuanSunCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelRumpfCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelSetschenow.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelEUNIQUAC.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelPitzerHMW.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct AqueousPhase::Impl
{
    /// The aqueous mixture instance
    AqueousMixture mixture;

    /// The base chemical model of the phase (yet to be combined with the custom activity coefficient models below)
    PhaseChemicalModel base_model;

    /// The functions that calculate the ln activity coefficients of selected species
    std::map<Index, AqueousActivityModel> ln_activity_coeff_functions;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance
    Impl(const AqueousMixture& mixture)
    : mixture(mixture)
    {}

    /// Return the combined chemical model function of the phase
    auto combinedChemicalModel() -> PhaseChemicalModel
    {
        // Create a copy of the data member `mixture` to be used in the following lambda function
        auto mixture = this->mixture;

        // Create a copy of the data member `base_model` to be used in the following lambda function
        auto base_model = this->base_model;

        // Create a copy of the data member `ln_activity_coeff_functions` to be used in the following lambda function
        auto ln_activity_coeff_functions = this->ln_activity_coeff_functions;

        // Define the function that calculates the chemical properties of the phase
        PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
        {
            // Evaluate the state of the aqueous mixture
            auto state = mixture.state(T, P, n);

            // Evaluate the aqueous chemical model
            base_model(res, T, P, n);
            
            // Update the activity coefficients and activities of selected species
            for(auto pair : ln_activity_coeff_functions)
            {
                const Index& i = pair.first; // the index of the selected species
                const AqueousActivityModel& func = pair.second; // the ln activity coefficient function of the selected species
                const ChemicalScalar ln_gi = func(state); // evaluate the ln activity coefficient function
                const ChemicalScalar ln_mi = log(state.m[i]); // get the molality of the selected species
                res.ln_activity_coefficients[i] = ln_gi; // update the ln activity coefficient selected species
                res.ln_activities[i] = ln_gi + ln_mi; // update the ln activity of the selected species
            }
        };

        return model;
    }
};

AqueousPhase::AqueousPhase()
: Phase(), pimpl(new Impl())
{
    setName("Aqueous");
    setType(PhaseType::Liquid);
}

AqueousPhase::AqueousPhase(const AqueousMixture& mixture)
: pimpl(new Impl(mixture))
{
    // Convert the AqueousSpecies instances to Species instances
    std::vector<Species> species;
    for(const AqueousSpecies& x : mixture.species())
        species.push_back(x);

    // Set the Phase attributes
    setName("Aqueous");
    setType(PhaseType::Liquid);
    setSpecies(species);
    setChemicalModelHKF();
    setActivityModelDuanSunCO2();
}

auto AqueousPhase::setInterpolationPoints(const std::vector<double>& temperatures, const std::vector<double>& pressures) -> AqueousPhase&
{
    pimpl->mixture.setInterpolationPoints(temperatures, pressures);
    return *this;
}

auto AqueousPhase::setChemicalModelIdeal() -> AqueousPhase&
{
    pimpl->ln_activity_coeff_functions.clear();
    pimpl->base_model = aqueousChemicalModelIdeal(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setChemicalModelDebyeHuckel() -> AqueousPhase&
{
	return setChemicalModelDebyeHuckel({});
}

auto AqueousPhase::setChemicalModelDebyeHuckel(const DebyeHuckelParams& params) -> AqueousPhase&
{
    pimpl->ln_activity_coeff_functions.clear();
    pimpl->base_model = aqueousChemicalModelDebyeHuckel(mixture(), params);
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setChemicalModelEUNIQUAC() -> AqueousPhase&
{
    return setChemicalModelEUNIQUAC({});
}

auto AqueousPhase::setChemicalModelEUNIQUAC(const EUNIQUACParams& params) -> AqueousPhase&
{
    pimpl->ln_activity_coeff_functions.clear();
    pimpl->base_model = aqueousChemicalModelEUNIQUAC(mixture(), params);
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setChemicalModelHKF() -> AqueousPhase&
{
    pimpl->ln_activity_coeff_functions.clear();
    pimpl->base_model = aqueousChemicalModelHKF(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setChemicalModelPitzerHMW() -> AqueousPhase&
{
    pimpl->ln_activity_coeff_functions.clear();
    pimpl->base_model = aqueousChemicalModelPitzerHMW(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModel(std::string species, const AqueousActivityModel& activity) -> AqueousPhase&
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = activity;
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModelIdeal(std::string species) -> AqueousPhase&
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelSetschenow(mixture(), 0.0);
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModelSetschenow(std::string species, double b) -> AqueousPhase&
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelSetschenow(mixture(), b);
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModelDuanSunCO2() -> AqueousPhase&
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelDuanSunCO2(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModelDrummondCO2() -> AqueousPhase&
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelDrummondCO2(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::setActivityModelRumpfCO2() -> AqueousPhase&
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelRumpfCO2(mixture());
    setChemicalModel(pimpl->combinedChemicalModel());
    return *this;
}

auto AqueousPhase::mixture() const -> const AqueousMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
