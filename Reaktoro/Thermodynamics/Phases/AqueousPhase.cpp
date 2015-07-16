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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Models/AqueousChemicalModelPitzerHMW.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityDrummond.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityDuanSun.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityIdeal.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityRumpf.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivitySetschenow.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct AqueousPhase::Impl
{
    /// The aqueous mixture instance
    AqueousMixture mixture;

    /// The functions that calculates the activities of selected species
    std::map<Index, AqueousActivityFunction> activities;

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
        // Define the function that calculates the chemical properties of the phase
        PhaseChemicalModel model = [=](double T, double P, const Vector& n)
        {
            // Calculate the state of the aqueous mixture
            const AqueousMixtureState state = mixture.state(T, P, n);

            // Evaluate the aqueous chemical model
            PhaseChemicalModelResult res = original(T, P, n);

            // Update the activity coefficients and activities of selected species
            for(auto pair : activities)
            {
                const Index& i = pair.first; // the index of the selected species
                const AqueousActivityFunction& ln_gi_func = pair.second; // the activity coefficient function of the selected species
                const ChemicalScalar ln_gi = ln_gi_func(state); // evaluate the activity coefficient function
                const ChemicalScalar ln_mi = log(state.m[i]); // get the molality of the selected species
                res.ln_activity_coefficients[i] = ln_gi; // update the activity coefficient selected species
                res.ln_activities[i] = ln_gi + ln_mi; // update the activity of the selected species
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

auto AqueousPhase::setActivityModel(std::string species, const AqueousActivityFunction& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = activity;
}

auto AqueousPhase::setActivityModelIdeal(std::string species) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = aqueousActivityIdeal(species, mixture());
}

auto AqueousPhase::setActivityModelSetschenow(std::string species, double b) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = aqueousActivitySetschenow(species, mixture(), b);
}

auto AqueousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = aqueousActivityDuanSunCO2(mixture());
}

auto AqueousPhase::setActivityModelDrummondCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = aqueousActivityDrummondCO2(mixture());
}

auto AqueousPhase::setActivityModelRumpfCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = aqueousActivityRumpfCO2(mixture());
}

auto AqueousPhase::mixture() const -> const AqueousMixture&
{
    return pimpl->mixture;
}

// todo delete these comments
//auto AqueousPhase::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
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
//auto AqueousPhase::activityConstants(double T, double P) const -> ThermoVector
//{
//    ThermoVector res(numSpecies());
//    res.val.setConstant(1.0);
//    return res;
//}
//
//auto AqueousPhase::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
//{
//    return activities(T, P, n)/concentrations(T, P, n);
//}
//
//auto AqueousPhase::activities(double T, double P, const Vector& n) const -> ChemicalVector
//{
//    AqueousMixtureState mixture_state = state(T, P, n);
//    const unsigned nspecies = numSpecies();
//    ChemicalVector a(nspecies, nspecies);
//    for(unsigned i = 0; i < nspecies; ++i)
//        a.row(i) = activity_fns[i](mixture_state);
//    return a;
//}

} // namespace Reaktoro
