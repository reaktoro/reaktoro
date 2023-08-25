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

#include "ActivityModelPhreeqc.hpp"

// C++ includes
#include <cassert>
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcUtils.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::log10;
using std::exp;

auto const LOG_10 = ln10;

auto createActivityModelPhreeqc(SpeciesList const& species, PhreeqcDatabase const& db) -> ActivityModel
{
    // Create the aqueous solution
    AqueousMixture solution(species);

    // The index of water in the solution
    auto const iw = solution.indexWater();

    // The molar mass of water (in kg/mol)
    auto const Mw = solution.species(iw).molarMass();

    // The PHREEQC pointer that gives us access to the internal state of PHREEQC after parsing the database.
    auto const phq = db.ptr();

    // The PHREEQC pointer to H2O species
    PhreeqcSpecies const* s_h2o = phq->s_h2o;

    // The llnl arrays in PHREEQC. These arrays are used to compute the activity coefficients of the species when using LLNL activity model.
    Vec<double> llnl_temp(phq->llnl_temp, phq->llnl_temp + phq->llnl_count_temp);
    Vec<double> llnl_adh(phq->llnl_adh, phq->llnl_adh + phq->llnl_count_adh);
    Vec<double> llnl_bdh(phq->llnl_bdh, phq->llnl_bdh + phq->llnl_count_bdh);
    Vec<double> llnl_bdot(phq->llnl_bdot, phq->llnl_bdot + phq->llnl_count_bdot);
    Vec<double> llnl_co2_coefs(phq->llnl_co2_coefs, phq->llnl_co2_coefs + phq->llnl_count_co2_coefs);

    // The PHREEQC species pointers for all species in the aqueous solution.
    Vec<PhreeqcSpecies const*> s_x;
    for(auto const& species : solution.species())
    {
        PhreeqcSpecies const* s = PhreeqcUtils::findSpecies(*phq, species.name());
        errorif(s == nullptr, "The species " + species.name() + " was not found in the PHREEQC database. Make sure you are assigning ActivityModelPhreeqc to an aqueous phase whose all species can be found in your chosen PHREEQC database.");
        s_x.push_back(s);
    }

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto aqstateptr = std::make_shared<AqueousMixtureState>();
    auto aqsolutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        // Check the validity of the mole fractions of the species (only at debug mode!)
        assert(x.minCoeff() > 0.0 && x.maxCoeff() <= 1.0);

        // Evaluate the state of the aqueous solution
        auto const& aqstate = *aqstateptr = solution.state(T, P, x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // Export the aqueous solution and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = aqstateptr;
        props.extra["AqueousMixture"] = aqsolutionptr;

        // Calculates gammas and [moles * d(ln gamma)/d mu] for all aqueous species.
        int i, j;
        int ifirst, ilast;
        real f, log_g_co2, dln_g_co2, c2_llnl;

        real c1, c2, a, b;
        real muhalf, equiv;
        real a_llnl, b_llnl, bdot_llnl;

        real mu = aqstate.Ie;
        real tk_x = T; // temperature in K
        real tc_x = T - 273.15; // temperature from K to °C
        real patm_x = P / 101325.0; // pressure from Pa to atm

        ArrayXr lg(s_x.size());

        // Initialize
        if(mu <= 0)
            mu = 1e-10;

        a_llnl = b_llnl = bdot_llnl = log_g_co2 = dln_g_co2 = c2_llnl = 0;

        // Compute temperature dependence of a and b for debye-huckel
        PhreeqcUtils::PhreeqcWaterProps wprops = PhreeqcUtils::waterPropsMemoized(T, P);
        a = wprops.wep.DH_A;
        b = wprops.wep.DH_B;

        // LLNL temperature dependence

        if(llnl_temp.size() > 0)
        {
            ifirst = 0;
            ilast = (int)llnl_temp.size();

            if(tc_x < llnl_temp[0] || tc_x > llnl_temp[llnl_temp.size() - 1])
                errorif(true, "Temperature out of range of LLNL_AQUEOUS_MODEL parameters");

            for(i = 0; i < (int)llnl_temp.size(); i++)
            {
                if(tc_x >= llnl_temp[i])
                    ifirst = i;
                if(tc_x <= llnl_temp[i])
                {
                    ilast = i;
                    break;
                }
            }
            if(ilast == ifirst)
            {
                f = 1;
            }
            else
            {
                f = (tc_x - llnl_temp[ifirst]) / (llnl_temp[ilast] - llnl_temp[ifirst]);
            }
            a_llnl = (1 - f) * llnl_adh[ifirst] + f * llnl_adh[ilast];
            b_llnl = (1 - f) * llnl_bdh[ifirst] + f * llnl_bdh[ilast];
            bdot_llnl = (1 - f) * llnl_bdot[ifirst] + f * llnl_bdot[ilast];

            // CO2 activity coefficient
            log_g_co2 = (llnl_co2_coefs[0] + llnl_co2_coefs[1] * tk_x + llnl_co2_coefs[2] / tk_x) * mu - (llnl_co2_coefs[3] + llnl_co2_coefs[4] * tk_x) * (mu / (mu + 1));
            log_g_co2 /= LOG_10;
            dln_g_co2 = (llnl_co2_coefs[0] + llnl_co2_coefs[1] * tk_x + llnl_co2_coefs[2] / tk_x) - (llnl_co2_coefs[3] + llnl_co2_coefs[4] * tk_x) * (1 / ((mu + 1) * (mu + 1)));
        }

        // constants for equations

        muhalf = sqrt(mu);
        c1 = (-a) * LOG_10 * (1.0 / (2 * muhalf * (muhalf + 1.0) * (muhalf + 1.0)) - 0.3);
        c2 = -a / (2 * muhalf);

        if(llnl_temp.size() > 0)
        {
            c2_llnl = -a_llnl / (2 * muhalf);
        }

        // Calculate activity coefficients
        for(i = 0; i < (int)s_x.size(); i++)
        {
            switch(s_x[i]->gflag)
            {
            // uncharged
            case 0:
                lg[i] = s_x[i]->dhb * mu;
                break;
            // Davies
            case 1:
                lg[i] = -s_x[i]->z * s_x[i]->z * a * (muhalf / (1.0 + muhalf) - 0.3 * mu);
                break;
            // Extended D-H, WATEQ D-H
            case 2:
                lg[i] = -a * muhalf * s_x[i]->z * s_x[i]->z / (1.0 + s_x[i]->dha * b * muhalf) + s_x[i]->dhb * mu;
                break;
            // Always 1.0
            case 3:
                lg[i] = 0.0;
                break;
            // Exchange
            case 4:
                errorif(true, "Exchange species should not exist in the aqueous phase.");
                break;
            // Always 1.0
            case 5:
                lg[i] = 0.0;
                break;
            // Surface
            case 6:
                errorif(true, "Surface species should not exist in the aqueous phase.");
                break;
            // LLNL
            case 7:
                if(llnl_temp.size() > 0)
                {
                    if(s_x[i]->z == 0)
                    {
                        lg[i] = 0.0;
                    }
                    else
                    {
                        lg[i] = -a_llnl * muhalf * s_x[i]->z * s_x[i]->z / (1.0 + s_x[i]->dha * b_llnl * muhalf) + bdot_llnl * mu;
                        break;
                    }
                }
                else
                {
                    errorif(true, "LLNL_AQUEOUS_MODEL_PARAMETERS not defined.");
                }
                break;
            // LLNL CO2
            case 8:
                if(llnl_temp.size() > 0)
                {
                    lg[i] = log_g_co2;
                }
                else
                {
                    errorif(true, "LLNL_AQUEOUS_MODEL_PARAMETERS not defined.");
                }
                break;
            // Activity coefficient of water
            case 9:
                lg[i] = 0.0; // Computed below
                break;
            }
        }

        // The mole fraction of water
        auto const xw = x[iw];

        // Set the activity coefficients of the species
        props.ln_g = lg * ln10;

        // Set the activities of the species
        props.ln_a = props.ln_g + log(aqstate.m);

        // Set the activitiy of water
        props.ln_a[iw] = log(1.0 - 0.017/Mw * (1.0 - xw)/xw); // From Equation (18) in PHREEQC's v2 manual; Note that: 1 - 0.0017*sum(n[i]/Waq) ≡ 1 - 0.017/Mw*sum(n[i]/n[w]) ≡ 1 - 0.017/Mw*sum(x[i]/x[w]) ≡ 1 - 0.017/Mw*(1 - x[w])/x[w] where sum is over solutes only; no water

        // Set the activity coefficient of water (mole fraction scale)
        props.ln_g[iw] = props.ln_a[iw] - log(xw);
    };

    return fn;
}

auto ActivityModelPhreeqc(PhreeqcDatabase const& db) -> ActivityModelGenerator
{
    // The PHREEQC pointer that gives us access to the internal state of PHREEQC after parsing the database.
    auto const phq = db.ptr();

    // Check if database object has been properly initialized
    errorif(phq == nullptr, "The PhreeqcDatabase object you are using with ActivityModelPhreeqc does not seem to have been properly initialized. Ensure you have specified a PHREEQC database file for its initialization.");

    // Check the database has a valid species pointer for water
    errorif(phq->s_h2o == nullptr, "The species pointer s_h2o in PHREEQC does not seem to have been initialized properly. Ensure the PHREEQC database you are using is valid and it contains a species that PHREEQC can recognize as water.");

    // Check if either a SIT or Pitzer based database is being used
    errorif(phq->sit_model, "ActivityModelPhreeqc currently does not support the SIT activity model of PHREEQC. Get in touch if this is a feature you need.");
    errorif(phq->pitzer_model, "ActivityModelPhreeqc currently does not support the Pitzer activity model of PHREEQC. Please use ActivityModelPitzer.");

    return [=](SpeciesList const& species) { return createActivityModelPhreeqc(species, db); };
}

} // namespace Reaktoro
