// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "AqueousChemicalModelHKF.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

/// The effective electrostatic radii of ionic species (in units of angstrom).
/// This data was taken from Table 3 of Helgeson et al. (1981).
const std::map<std::string, double> effective_radii =
{
    {"H+"  , 3.08}, {"Fe+++", 3.46},
    {"Li+" , 1.64}, {"Al+++", 3.33},
    {"Na+" , 1.91}, {"Au+++", 3.72},
    {"K+"  , 2.27}, {"La+++", 3.96},
    {"Rb+" , 2.41}, {"Gd+++", 3.79},
    {"Cs+" , 2.61}, {"In+++", 3.63},
    {"NH4+", 2.31}, {"Ca+++", 3.44},
    {"Ag+" , 2.20}, {"F-"   , 1.33},
    {"Au+" , 2.31}, {"Cl-"  , 1.81},
    {"Cu+" , 1.90}, {"Br-"  , 1.96},
    {"Mg++", 2.54}, {"I-"   , 2.20},
    {"Sr++", 3.00}, {"OH-"  , 1.40},
    {"Ca++", 2.87}, {"HS-"  , 1.84},
    {"Ba++", 3.22}, {"NO3-" , 2.81},
    {"Pb++", 3.08}, {"HCO3-", 2.10},
    {"Zn++", 2.62}, {"HSO4-", 2.37},
    {"Cu++", 2.60}, {"ClO4-", 3.59},
    {"Cd++", 2.85}, {"ReO4-", 4.23},
    {"Hg++", 2.98}, {"SO4--", 3.15},
    {"Fe++", 2.62}, {"CO3--", 2.81},
    {"Mn++", 2.68}
};

/// Calculate the effective electrostatic radius of species (in units of A).
/// @param species The aqueous species instance of the ionic species
/// @return The effective electrostatic radius of the ionic species (in units of A)
double effectiveIonicRadius(const AqueousSpecies& species)
{
    // Find the effective ionic radius of the species in `effective_radii`.
    // Note that `species` might have a different name convention than those
    // used in `effective_radii`. Thus, we need to check if the name of given
    // `species` is an alternative to a name in `effective_radii`.
    for(auto pair : effective_radii)
        if(isAlternativeChargedSpeciesName(species.name(), pair.first))
            return pair.second;

    // The electrical charge of the species
    const double Zi = species.charge();

    // Estimated effective ionci radius of the species based on TOUGHREACT approach
    if(Zi == -1) return 1.81;        // based on Cl- value
    if(Zi == -2) return 3.00;        // based on rounded average of CO3-- and SO4-- values
    if(Zi == -3) return 4.20;        // based on estimation from straight line fit with charge
    if(Zi == +1) return 2.31;        // based on NH4+ value
    if(Zi == +2) return 2.80;        // based on rounded average for +2 species in the HKF table of effective ionic radii
    if(Zi == +3) return 3.60;        // based on rounded average for +3 species in the HKF table of effective ionic radii
    if(Zi == +4) return 4.50;        // based on estimaton using HKF eq. 142
    if(Zi <  -3) return -Zi*4.2/3.0; // based on linear extrapolation
    return Zi*4.5/4.0;               // based on linear extrapolation
}

} // namespace
auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture) -> PhaseChemicalModel
{
    // The number of species in the mixture
    const unsigned num_species = mixture.numSpecies();

    // The number of charged species in the mixture
    const unsigned num_charged_species = mixture.numChargedSpecies();

    // The indices of the charged species
    const Indices icharged_species = mixture.indicesChargedSpecies();

    // The index of the water species
    const Index iwater = mixture.indexWater();

    // The effective electrostatic radii of the charged species
    std::vector<double> effective_radii;

    // The electrical charges of the charged species only
    std::vector<double> charges;

    // The natural log of 10
    const double ln10 = std::log(10);

    // The molar mass of water
    const double Mw = waterMolarMass;

    // Collect the effective radii of the ions
    for(Index idx_ion : icharged_species)
    {
        const AqueousSpecies& species = mixture.species(idx_ion);
        effective_radii.push_back(effectiveIonicRadius(species));
        charges.push_back(species.charge());
    }

    // Define the intermediate chemical model function of the aqueous mixture
    auto model = [=](const AqueousMixtureState& state)
    {
        // Auxiliary references to state variables
        const auto& T = state.T;
        const auto& I = state.Ie;
        const auto& x = state.x;
        const auto& m = state.m;
        const auto& rho = state.rho/1000; // density in g/cm3
        const auto& epsilon = state.epsilon;

        // The molar fraction of the water species and its molar derivatives
        const auto xw = x.row(iwater);

        // The ln and log10 of water molar fraction
        const auto ln_xw = log(xw);
        const auto log10_xw = log10(xw);

        // The square root of the ionic strength
        const auto sqrtI = sqrt(I);
        const auto sqrt_rho = sqrt(rho);
        const auto sqrt_T_epsilon = sqrt(T * epsilon);
        const auto T_epsilon = T * epsilon;

        // The parameters for the HKF model
        const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
        const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;

        // The alpha parameter used in the calculation of osmotic coefficient of water
        const auto alpha = xw/(1.0 - xw) * log10_xw;

        // The osmotic coefficient of the aqueous phase
        ChemicalScalar phi(num_species);

        // The result of the equation of state
        PhaseChemicalModelResult res(num_species);

        // Set the activity coefficients of the neutral species to
        // water molar fraction to convert it to molality scale
        res.ln_activity_coefficients = ln_xw;

        // Loop over all charged species in the mixture
        for(unsigned i = 0; i < num_charged_species; ++i)
        {
            // The index of the charged species in the mixture
            const Index ispecies = icharged_species[i];

            // The molality of the charged species and its molar derivatives
            const auto mi = m.row(ispecies);

            // Check if the molality of the charged species is zero
            if(mi.val == 0.0)
                continue;

            // The electrical charge of the charged species
            const auto z = charges[i];
            const auto z2 = z*z;

            // The effective radius of the charged species
            const auto eff_radius = effective_radii[i];

            // The Debye-Huckel ion size parameter of the current ion as computed by Reed (1982) and also in TOUGHREACT
            const auto a = (z < 0) ?
                2.0*(eff_radius + 1.91*std::abs(z))/(std::abs(z) + 1.0) :
                2.0*(eff_radius + 1.81*std::abs(z))/(std::abs(z) + 1.0);

            // The \Lamba parameter of the HKF activity coefficient model and its molar derivatives
            const ChemicalScalar lambda = 1.0 + a*B*sqrtI;

            // The log10 activity coefficient of the charged species (in molality scale) and its molar derivatives
            const auto log10_gi = -(A*z2*sqrtI)/lambda + log10_xw;

            // Set the activity coefficient of the current charged species
            res.ln_activity_coefficients[ispecies] = log10_gi * ln10;

            // Check if the molar fraction of water is one
            if(xw != 1.0)
            {
                // The sigma parameter of the current ion and its molar derivatives
                const auto sigma = 3.0/pow(a*B*sqrtI, 3) * (lambda - 1.0/lambda - 2.0*log(lambda));

                // The psi contribution of the current ion and its molar derivatives
                const auto psi = A*z2*sqrtI*sigma/3.0 + alpha;

                // Update the osmotic coefficient with the contribution of the current charged species
                phi += mi * psi;
            }
        }

        // Set the activities of the solutes (molality scale)
        res.ln_activities = res.ln_activity_coefficients + log(m);

        // Set the activity of water (in molar fraction scale)
        if(xw != 1.0) res.ln_activities[iwater] = ln10 * Mw * phi;
                 else res.ln_activities[iwater] = ln_xw;

        // Set the activity coefficient of water (molar fraction scale)
        res.ln_activity_coefficients[iwater] = res.ln_activities[iwater] - ln_xw;

        // Set the activity constants of aqueous species to ln(55.508472)
        res.ln_activity_constants = std::log(55.508472);

        // Set the activity constant of water to zero
        res.ln_activity_constants[iwater] = 0.0;

        return res;
    };

    // Define the chemical model function of the aqueous mixture
    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        // Calculate the state of the mixture
        const AqueousMixtureState state = mixture.state(T, P, n);

        return model(state);
    };

    return f;
}

DebyeHuckelParams::DebyeHuckelParams()
{}


} // namespace Reaktoro
