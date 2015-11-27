// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "AqueousPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDrummondCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDuanSunCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelRumpfCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelSetschenow.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelPitzerHMW.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct AqueousPhase::Impl
{
    /// The aqueous mixture instance
    AqueousMixture mixture;

    /// The functions that calculate the ln activity coefficients of selected species
    std::map<Index, AqueousActivityModel> ln_activity_coeff_functions;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance
    Impl(const AqueousMixture& mixture)
    : mixture(mixture)
    {}

    /// Return the chemical model function of the phase
    auto convertPhaseChemicalModel(const PhaseChemicalModel& original) -> PhaseChemicalModel
    {
        // Create a copy of the data member `mixture` to be used in the following lambda function
        auto mixture = this->mixture;

        // Create a copy of the data member `ln_activity_coeff_functions` to be used in the following lambda function
        auto ln_activity_coeff_functions = this->ln_activity_coeff_functions;

        // Define the function that calculates the chemical properties of the phase
        PhaseChemicalModel model = [=](Temperature T, Pressure P, const Vector& n)
        {
            // Calculate the state of the aqueous mixture
            const AqueousMixtureState state = mixture.state(T, P, n);

            // Evaluate the aqueous chemical model
            PhaseChemicalModelResult res = original(T, P, n);

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

            return res;
        };

        return model;
    }
};

AqueousPhase::AqueousPhase()
: pimpl(new Impl())
{}

AqueousPhase::AqueousPhase(const AqueousPhase& other)
: Phase(other), pimpl(new Impl(*other.pimpl))
{}

AqueousPhase::AqueousPhase(const AqueousMixture& mixture)
: pimpl(new Impl(mixture))
{
    // Convert the AqueousSpecies instances to Species instances
    std::vector<Species> species;
    for(const AqueousSpecies& x : mixture.species())
        species.push_back(x);

    // Set the Phase attributes
    setName("Aqueous");
    setSpecies(species);
    setReferenceState(PhaseReferenceState::IdealSolution);
    setChemicalModelHKF();
    setActivityModelDuanSunCO2();
}

AqueousPhase::~AqueousPhase()
{}

auto AqueousPhase::operator=(AqueousPhase other) -> AqueousPhase&
{
    Phase::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousPhase::setChemicalModelIdeal() -> void
{
    // Create the aqueous chemical model
    PhaseChemicalModel aqueous_model = aqueousChemicalModelIdeal(mixture());

    // Convert the PhaseChemicalModel to PhaseChemicalModel
    PhaseChemicalModel model = pimpl->convertPhaseChemicalModel(aqueous_model);

    setChemicalModel(model);
}

auto AqueousPhase::setChemicalModelHKF() -> void
{
    // Create the aqueous chemical model
    PhaseChemicalModel aqueous_model = aqueousChemicalModelHKF(mixture());

    // Convert the PhaseChemicalModel to PhaseChemicalModel
    PhaseChemicalModel model = pimpl->convertPhaseChemicalModel(aqueous_model);

    setChemicalModel(model);
}

auto AqueousPhase::setChemicalModelPitzerHMW() -> void
{
    // Create the aqueous chemical model
    PhaseChemicalModel aqueous_model = aqueousChemicalModelPitzerHMW(mixture());

    // Convert the PhaseChemicalModel to PhaseChemicalModel
    PhaseChemicalModel model = pimpl->convertPhaseChemicalModel(aqueous_model);

    setChemicalModel(model);
}

auto AqueousPhase::setActivityModel(std::string species, const AqueousActivityModel& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = activity;
}

auto AqueousPhase::setActivityModelIdeal(std::string species) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelSetschenow(mixture(), 0.0);
}

auto AqueousPhase::setActivityModelSetschenow(std::string species, double b) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelSetschenow(mixture(), b);
}

auto AqueousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelDuanSunCO2(mixture());
}

auto AqueousPhase::setActivityModelDrummondCO2() -> void
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelDrummondCO2(mixture());
}

auto AqueousPhase::setActivityModelRumpfCO2() -> void
{
    const Index ispecies = indexSpeciesAny(alternativeNeutralSpeciesNames("CO2(aq)"));
    if(ispecies < numSpecies())
        pimpl->ln_activity_coeff_functions[ispecies] = aqueousActivityModelRumpfCO2(mixture());
}

auto AqueousPhase::mixture() const -> const AqueousMixture&
{
    return pimpl->mixture;
}

// todo delete these comments
//auto AqueousPhase::concentrations(Temperature T, Pressure P, const Vector& n) const -> ChemicalVector
//{
//    // Calculate the molalities of the species
//    ChemicalVector c = molalities(n);
//
//    // Calculate the molar fractions of the species
//    ChemicalVector x = molarFractions(n);
//
//    // The index of the water species
//    const Index iH2O = indexWater();
//
//    // Set the concentration of water to its molar fraction
//    c.row(iH2O) = x.row(iH2O);
//
//    return c;
//}
//
//auto AqueousPhase::activityConstants(Temperature T, Pressure P) const -> ThermoVector
//{
//    ThermoVector res(numSpecies());
//    res.val.setConstant(1.0);
//    return res;
//}
//
//auto AqueousPhase::activityCoefficients(Temperature T, Pressure P, const Vector& n) const -> ChemicalVector
//{
//    return activities(T, P, n)/concentrations(T, P, n);
//}
//
//auto AqueousPhase::activities(Temperature T, Pressure P, const Vector& n) const -> ChemicalVector
//{
//    AqueousMixtureState mixture_state = state(T, P, n);
//    const unsigned nspecies = numSpecies();
//    ChemicalVector a(nspecies, nspecies);
//    for(unsigned i = 0; i < nspecies; ++i)
//        a.row(i) = activity_fns[i](mixture_state);
//    return a;
//}

} // namespace Reaktoro
