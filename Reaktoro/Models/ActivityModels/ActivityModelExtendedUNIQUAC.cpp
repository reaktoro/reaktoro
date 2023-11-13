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

#include "ActivityModelExtendedUNIQUAC.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>

namespace Reaktoro {

auto createActivityModelExtendedUNIQUAC(SpeciesList const& species, ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModel
{
    // Create the aqueous solution and its state object
    AqueousMixture solution(species);
    AqueousMixtureState aqstate;

    // The electrical charges of the species
    const ArrayXd z = solution.charges();

    // The electrical charges of the species squared
    const ArrayXd z2 = z * z;

    // The index of water in the solution
    const Index iH2O = solution.indexWater();

    Indices iDH; // The indices of the species containing Debye-Huckel contribution
    Indices iR; // The indices of the species containing the surface area parameter rᵢ
    Indices iQ; // The indices of the species containing the volume parameter qᵢ
    Indices iU; // The indices of the species containing energetic binary interaction parameters uᵢⱼ
    // Indices iUQ; // The indices of the species containing both the volume parameter qᵢ and the energetic binary interaction parameters uᵢⱼ

    ArrayXr r = ArrayXr::Zero(species.size()); // TODO: Initialize properly with values from parameters
    ArrayXr q = ArrayXr::Zero(species.size()); // TODO: Initialize properly with values from parameters

    MatrixXr u0 = MatrixXr::Constant(species.size(), species.size(), 1.0e+16); // The matrix u with the energy interaction parameters in the e-UNIQUAC model
    MatrixXr uT = MatrixXr::Constant(species.size(), species.size(), 1.0e+16); // The matrix u with the energy interaction parameters in the e-UNIQUAC model

    u0.diagonal().fill(0.0);
    uT.diagonal().fill(0.0);

    ArrayXr lnGammaDH(species.size());
    ArrayXr lnGammaCombinatorial(species.size());
    ArrayXr lnGammaResidual(species.size());

    ArrayXr thetapsi;  // The matrix-vector product tr(ψ)θ

    for(auto const& [formula, param] : params.r)
    {
        errorif(param <= 0.0, "The surface area parameter rᵢ in the Extended UNIQUAC model cannot be zero or negative: r[", formula, "] = ", param);
        const auto ispecies = species.findWithFormula(formula);
        if(ispecies >= species.size())
            continue;
        iR.push_back(ispecies);
        r[ispecies] = param;
    }

    for(auto const& [formula, param] : params.q)
    {
        errorif(param <= 0.0, "The volume parameter qᵢ in the Extended UNIQUAC model cannot be zero or negative: q[", formula, "] = ", param);
        const auto ispecies = species.findWithFormula(formula);
        if(ispecies >= species.size())
            continue;
        iQ.push_back(ispecies);
        q[ispecies] = param;
    }

    for(auto const& [formula1, formula2, params] : params.u)
    {
        const auto ispecies1 = species.findWithFormula(formula1);
        const auto ispecies2 = species.findWithFormula(formula2);
        if(ispecies1 >= species.size() || ispecies2 >= species.size())
            continue;
        iU.push_back(ispecies1);
        iU.push_back(ispecies2);
        u0(ispecies2, ispecies1) = params[0];
        uT(ispecies2, ispecies1) = params[1];
    }

    iR = unique(iR);
    iQ = unique(iQ);
    iU = unique(iU);
    // iUQ = intersect(iU, iQ);

    // errorif(identical(iR, iQ), "Expecting both rᵢ and qᵢ parameters in the Extended UNIQUAC model, but for some species, only one of these parameters were provided.");

    r = ArrayXr(r(iR));
    q = ArrayXr(q(iQ));

    u0 = MatrixXr(u0(iQ, iQ)); // Only binary interaction parameters involving species with positive qᵢ values
    uT = MatrixXr(uT(iQ, iQ)); // Only binary interaction parameters involving species with positive qᵢ values

    ArrayXr phi; // The ϕ array in the e-UNIQUAC model
    ArrayXr theta; // The θ array in the e-UNIQUAC model

    MatrixXr u = MatrixXr::Zero(iQ.size(), iQ.size()); // The matrix u with the energy interaction parameters in the e-UNIQUAC model
    MatrixXr psi = MatrixXr::Zero(iQ.size(), iQ.size()); // The matrix ψ in the e-UNIQUAC model

    ArrayXr xr;
    ArrayXr xq;

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& T = args.T;
        auto const& P = args.P;
        auto const& x = args.x;

        // Evaluate the state of the aqueous solution
        auto const& aqstate = solution.state(T, P, x);

        auto const& I = aqstate.Ie;
        auto const& sqrtI = sqrt(I);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // The mole fraction of water
        auto const xw = x[iH2O];

        // Calculate the Debye-Huckel A coefficient according to equation (6) of Thomsen (2005)
        const auto Tr = 273.15;
        const auto A = 1.131 + 1.335e-3*(T - Tr) + 1.164e-5*(T - Tr)*(T - Tr);
        const auto b = 1.5; // This is the aB term in the Debye-Huckel contribution (assumed constant in the e-UNIQUAC model)
        const auto Z = 10.0; // The co-ordination number z of the UNIQUAC model (Thomasen, 1997, equation 4.2)
        const auto Lambda = 1.0 + b*sqrtI;
        auto const alpha = A*sqrtI/Lambda;

        xr = x(iR);
        xq = x(iQ);

        phi = (xr * r) / sum(xr * r);
        theta = (xq * q) / sum(xq * q);

        for(auto j = 0; j < iQ.size(); ++j)
            for(auto i = j; i < iQ.size(); ++i)
                u(j, i) = u(i, j) = u0(j, i) + uT(j, i)*(T - Tr);

        for(auto j = 0; j < iQ.size(); ++j)
            for(auto i = 0; i < iQ.size(); ++i)
                psi(j, i) = exp(-(u(j, i) - u(i, i))/T);

        thetapsi = tr(psi) * theta.matrix();

        lnGammaDH(iDH) = -z2(iDH)*alpha;

        for(auto const& [i, ispecies] : enumerate(iR))
            lnGammaCombinatorial[ispecies] = log(phi[i]/x[i]) + 1 - phi[i]/x[i] - 0.5*Z*q[i]*(log(phi[i]/theta[i]) + 1 - phi[i]/theta[i]);

        for(auto const& [i, ispecies] : enumerate(iQ))
            lnGammaResidual[ispecies] = q[i]*(1 - log(thetapsi[i]) - (theta * ArrayXr(psi.row(i))/thetapsi).sum());

        // Set the activity coefficients of the solutes
        // props.ln_g = pzstate.ln_gamma;

        // // Set the activities of the solutes
        // props.ln_a = props.ln_g + log(aqstate.m);

        // // Set the activitiy of water
        // props.ln_a[iH2O] = pzstate.ln_aw;

        // // Set the activity coefficient of water (mole fraction scale)
        // props.ln_g[iH2O] = pzstate.ln_aw - log(xw);
    };

    return fn;
}

auto ActivityModelExtendedUNIQUAC() -> ActivityModelGenerator
{
    return ActivityModelExtendedUNIQUAC(Params::embedded("ExtendedUNIQUAC/Thomsen1997.yaml"));
}

auto ActivityModelExtendedUNIQUAC(ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelExtendedUNIQUAC(species, params); };
}

auto ActivityModelExtendedUNIQUAC(Params const& params) -> ActivityModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ActivityModelParams"), "Expecting ExtendedUNIQUAC activity model parameters in given Params object, but it lacks a `ActivityModelParams` section within which another section `ExtendedUNIQUAC` should exist.");
    errorif(!data.at("ActivityModelParams").exists("ExtendedUNIQUAC"), "Expecting ExtendedUNIQUAC activity model parameters in given Params object, under the section `ExtendedUNIQUAC`.");
    errorif(!data.at("ActivityModelParams").at("ExtendedUNIQUAC").isDict(), "Expecting section `ExtendedUNIQUAC` with ExtendedUNIQUAC activity model parameters to be a dictionary.");

    ActivityModelParamsExtendedUNIQUAC euniquacparams =
        data["ActivityModelParams"]["ExtendedUNIQUAC"].as<ActivityModelParamsExtendedUNIQUAC>();

    return ActivityModelExtendedUNIQUAC(euniquacparams);
}

} // namespace Reaktoro
