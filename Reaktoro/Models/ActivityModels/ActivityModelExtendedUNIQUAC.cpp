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

// C++ includes
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

using std::log;
using std::sqrt;
using std::exp;
using std::abs;
using std::round;

/// Return the EUNIQUAC temperature-pressure correction model based on constant expression.
auto createParamCorrectionModelConstant(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;

    errorif(c.size() == 0, "Cannot create the Constant temperature-dependent EUNIQUAC parameter function with empty coefficients");
    errorif(c.size() >= 1, "Cannot create the Constant temperature-dependent EUNIQUAC parameter function with more than 1 coefficient. Given coefficients: ", c);

    return [=](real const& T, real const& Pbar) { return c[0]; };
}

/// Return the EUNIQUAC temperature-pressure correction model based on a linear temperature function.
auto createParamCorrectionModelLinear(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;

    errorif(c.size() == 0, "Cannot create the Linear temperature-dependent EUNIQUAC parameter function with empty coefficients");
    errorif(c.size() >= 2, "Cannot create the Linear temperature-dependent EUNIQUAC parameter function with more than 2 coefficients. Given coefficients: ", c);

    if(c.size() == 1) return [=](real const& T, real const& Pbar) { return c[0]; };
    if(c.size() == 2) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*T; };
}

/// Return the EUNIQUAC temperature-pressure correction model based on a quadratic temperature function.
auto createParamCorrectionModelQuadratic(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;

    errorif(c.size() == 0, "Cannot create the Quadratic temperature-dependent EUNIQUAC parameter function with empty coefficients");
    errorif(c.size() >= 3, "Cannot create the Quadratic temperature-dependent EUNIQUAC parameter function with more than 3 coefficients. Given coefficients: ", c);

    if(c.size() == 1) return [=](real const& T, real const& Pbar) { return c[0]; };
    if(c.size() == 2) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*T; };
    if(c.size() == 3) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*T + c[2]*T*T; };
}

/// Return the EUNIQUAC temperature-pressure correction model with given coefficients and model option.
auto createParamCorrectionModel(Vec<Param> const& coefficients, ActivityModelParamsExtendedUNIQUAC::CorrectionModel const& option) -> Fn<real(real const&, real const&)>
{
    switch(option)
    {
        case ActivityModelParamsExtendedUNIQUAC::CorrectionModel::Constant: return createParamCorrectionModelConstant(coefficients);
        case ActivityModelParamsExtendedUNIQUAC::CorrectionModel::Linear: return createParamCorrectionModelLinear(coefficients);
        case ActivityModelParamsExtendedUNIQUAC::CorrectionModel::Quadratic: return createParamCorrectionModelQuadratic(coefficients);
        default: return createParamCorrectionModelConstant(coefficients);
    }
}

/// Used to represent an interaction parameter in the EUNIQUAC activity model.
struct PitzerParam
{
    /// The indices of the species associated with this interaction parameter.
    const Indices ispecies;

    /// The temperature-pressure correction model for the interaction parameter.
    const Fn<real(real const&, real const&)> model;

    /// The current value of the interaction parameter since last update.
    real value;

    /// Update the current value of the interaction parameter.
    auto update(real const& T, real const& P)
    {
        value = model(T, P);
    }
};

/// Auxiliary alias for ActivityModelParamsExtendedUNIQUAC::InteractionParamAttribs.
using InteractionParamAttribsEUNIQUAC = ActivityModelParamsExtendedUNIQUAC::InteractionParamAttribs;

/// Return a given list of species indices sorted in ascending order of electric charge.
auto sortedSpeciesIndicesByCharge(SpeciesList const& specieslist, Indices const& ispecies) -> Indices
{
    // Lambda function that returns the charge of the ith species in `specieslist`
    auto chargefn = [&](auto i) { return specieslist[i].charge(); };

    // Lambda function that returns true if two species indices are in ascending order of electric charge.
    auto lesschargefn = [&](Index ispecies1, Index ispecies2) { return chargefn(ispecies1) < chargefn(ispecies2); };

    return sortedfn(ispecies, lesschargefn);
}

/// Convert a InteractionParamAttribsEUNIQUAC object to a PitzerParam one for a binary interaction parameter.
auto createPitzerParamBinary(SpeciesList const& specieslist, InteractionParamAttribsEUNIQUAC const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 2, "Expecting two chemical formulas when processing a Pitzer binary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const numspecies = specieslist.size();

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    if(ispecies1 >= numspecies)
        return {}; // species1 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);
    if(ispecies2 >= numspecies)
        return {}; // species2 is not present in the list of aqueous species; return empty PitzerParam object

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.parameters, attribs.model) };
}

/// Convert a InteractionParamAttribsEUNIQUAC object to a PitzerParam one for a ternary interaction parameter.
auto createPitzerParamTernary(SpeciesList const& specieslist, InteractionParamAttribsEUNIQUAC const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 3, "Expecting two chemical formulas when processing a Pitzer ternary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const numspecies = specieslist.size();

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    if(ispecies1 >= numspecies)
        return {}; // species1 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);
    if(ispecies2 >= numspecies)
        return {}; // species2 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies3 = specieslist.findWithFormula(attribs.formulas[2]);
    if(ispecies3 >= numspecies)
        return {}; // species3 is not present in the list of aqueous species; return empty PitzerParam object

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2, ispecies3});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.parameters, attribs.model) };
}

/// The auxiliary type used to store computed values from the Pitzer activity model implemented by @ref Pitzer.
struct PitzerState
{
    real phiw;        ///< The osmotic coefficient of water.
    real ln_aw;       ///< The activity of water (natural log).
    ArrayXr ln_gamma; ///< The activity coefficients of the aqueous species (natural log).
};

/// The auxiliary type used to implement the Pitzer activity model.
struct PitzerModel
{
    AqueousMixture solution; ///< The aqueous solution for which this Pitzer activity model is defined.

    Vec<PitzerParam> beta0;  ///< The parameters \eq{\beta^{(0)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> beta1;  ///< The parameters \eq{\beta^{(1)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> beta2;  ///< The parameters \eq{\beta^{(2)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> Cphi;   ///< The parameters \eq{C^{\phi}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> theta;  ///< The parameters \eq{\theta_{ij}(T, P)} in the Pitzer model for cation-cation and anion-anion interactions.
    Vec<PitzerParam> psi;    ///< The parameters \eq{\psi_{ijk}(T, P)} in the Pitzer model for cation-cation-anion and anion-anion-cation interactions.
    Vec<PitzerParam> lambda; ///< The parameters \eq{\lambda_{ij}(T, P)} in the Pitzer model for neutral-cation and neutral-anion interactions.
    Vec<PitzerParam> zeta;   ///< The parameters \eq{\zeta_{ijk}(T, P)} in the Pitzer model for neutral-cation-anion interactions.
    Vec<PitzerParam> mu;     ///< The parameters \eq{\mu_{ijk}(T, P)} in the Pitzer model for neutral-neutral-neutral, neutral-neutral-cation, and neutral-neutral-anion interactions.
    Vec<PitzerParam> eta;    ///< The parameters \eq{\eta_{ijk}(T, P)} in the Pitzer model for neutral-cation-cation and neutral-anion-anion interactions.

    Vec<Param> alpha1; ///< The parameters \eq{alpha_1_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.
    Vec<Param> alpha2; ///< The parameters \eq{alpha_2_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.

    using Tuples2i = Tuples<Index, Index>;                   ///< Auxiliary type for a tuple of 2 index values.
    using Tuples3d = Tuples<double, double, double>;         ///< Auxiliary type for a tuple of 3 double values.
    using Tuples4d = Tuples<double, double, double, double>; ///< Auxiliary type for a tuple of 4 double values.

    Tuples3d lambda_coeffs; ///< The coefficients multiplying the terms where the lambda Pitzer parameter is involved.
    Tuples4d mu_coeffs;     ///< The coefficients multiplying the terms where the mu Pitzer parameter is involved.

    Tuples2i thetaij; ///< The indices (i, j) of the cation-cation and anion-anion species pairs used to account for \eq{^{E}\theta_{ij}(I)} and \eq{^{E}\theta_{ij}^{\prime}(I)} contributions.
    ArrayXr thetaE;   ///< The current values of the parameters \eq{^{E}\theta_{ij}(I)} associated to the \eq{\theta_{ij}} parameters.
    ArrayXr thetaEP;  ///< The current values of the parameters \eq{^{E}\theta_{ij}^{\prime}(I)} associated to the \eq{\theta_{ij}} parameters.

    Fn<real(real const&, real const&)> Aphi; ///< The function that computes the Debye-huckel parameter \eq{A^\phi(T, P)} in the Pitzer model.

    /// Construct a default Pitzer object.
    PitzerModel()
    {}

    /// Construct a Pitzer object with given list of species in the aqueous solution and the parameters for the Pitzer activity model.
    // PitzerModel(SpeciesList const& specieslist, ActivityModelParamsExtendedUNIQUAC const& params)
    PitzerModel(AqueousMixture const& solution, ActivityModelParamsExtendedUNIQUAC const& params)
    : solution(solution)
    {
        for(auto const& entry : params.beta0)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta0.push_back(param);

        for(auto const& entry : params.beta1)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta1.push_back(param);

        for(auto const& entry : params.beta2)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta2.push_back(param);

        for(auto const& entry : params.Cphi)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                Cphi.push_back(param);

        for(auto const& entry : params.theta)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                theta.push_back(param);

        for(auto const& entry : params.psi)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                psi.push_back(param);

        for(auto const& entry : params.lambda)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                lambda.push_back(param);

        for(auto const& entry : params.zeta)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                zeta.push_back(param);

        for(auto const& entry : params.mu)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                mu.push_back(param);

        for(auto const& entry : params.eta)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                eta.push_back(param);

        for(auto const& entry : params.beta1)
            alpha1.push_back(determineAlpha1(entry.formulas[0], entry.formulas[1], params.alpha1));

        for(auto const& entry : params.beta2)
            alpha2.push_back(determineAlpha2(entry.formulas[0], entry.formulas[1], params.alpha2));

        auto const& ications = solution.indicesCations();
        auto const& ianions = solution.indicesAnions();

        for(auto i = 0; i < ications.size() - 1; ++i)
            for(auto j = i + 1; j < ications.size(); ++j)
                thetaij.emplace_back(ications[i], ications[j]);

        for(auto i = 0; i < ianions.size() - 1; ++i)
            for(auto j = i + 1; j < ianions.size(); ++j)
                thetaij.emplace_back(ianions[i], ianions[j]);

        auto const& z = solution.charges();

        for(auto const& entry : lambda)
        {
            auto const i1 = entry.ispecies[0];
            auto const i2 = entry.ispecies[1];
            lambda_coeffs.push_back(determineLambdaCoeffs(z[i1], z[i2], i1, i2));
        }

        for(auto const& entry : mu)
        {
            auto const i1 = entry.ispecies[0];
            auto const i2 = entry.ispecies[1];
            auto const i3 = entry.ispecies[2];
            mu_coeffs.push_back(determineMuCoeffs(z[i1], z[i2], z[i3], i1, i2, i3));
        }

        // Define the function Aphi(T, P) according to PHREEQC (see method calc_dielectrics at utilities.cpp for computing A0)
        Aphi = [](real const& T, real const& P) -> real
        {
            auto const rho = PhreeqcUtils::waterDensityPhreeqc(T, P)/1000; // the density of water (in g/cm3)
            auto const epsilon = PhreeqcUtils::waterDielectricConstantPhreeqc(T, P);
            auto const AVOGADRO = 6.02252e23;
            auto const pi = 3.14159265358979;
            auto const e2_DkT = 1.671008e-3 / (epsilon * T);
            auto const DH_B = sqrt(8 * pi * AVOGADRO * e2_DkT * rho / 1e3); // Debye length parameter, 1/cm(mol/kg)^-0.5
            auto const Aphi0 = DH_B * e2_DkT / 6.0;
            return Aphi0;
        };

        Aphi = memoizeLast(Aphi); // memoize so that subsequent repeated calls with same (T, P) return cached result.
    }

    /// Update all Pitzer interaction parameters according to current temperature and pressure.
    auto updateParams(real const& T, real const& P)
    {
        // Note: Do not try to avoid the calculations before in case T and P is
        // the same as last time because T and P could be the same but the Param
        // objects in these models could be changing!

        auto const Pbar = P * 1e-5; // from Pa to bar

        for(auto& param : beta0)
            param.value = param.model(T, Pbar);

        for(auto& param : beta1)
            param.value = param.model(T, Pbar);

        for(auto& param : beta2)
            param.value = param.model(T, Pbar);

        for(auto& param : Cphi)
            param.value = param.model(T, Pbar);

        for(auto& param : theta)
            param.value = param.model(T, Pbar);

        for(auto& param : psi)
            param.value = param.model(T, Pbar);

        for(auto& param : lambda)
            param.value = param.model(T, Pbar);

        for(auto& param : zeta)
            param.value = param.model(T, Pbar);

        for(auto& param : mu)
            param.value = param.model(T, Pbar);

        for(auto& param : eta)
            param.value = param.model(T, Pbar);
    }

    /// Evaluate the Pitzer model and compute the properties of the aqueous solution.
    auto evaluate(AqueousMixtureState const& aqstate, PitzerState& pzstate)
    {
        auto const& T = aqstate.T; // in K
        auto const& P = aqstate.P; // in Pa
        auto const& M = aqstate.m; // in molal
        auto const& z = solution.charges();
        auto const& iH2O = solution.indexWater();
        auto const& icharged = solution.indicesCharged();
        auto const& Mw = solution.water().molarMass(); // in kg/mol

        /// Update all Pitzer interaction parameters according to current temperature and pressure.
        updateParams(T, P);

        // The ionic strength of the solution and its square-root
        auto const I = aqstate.Ie;
        auto const DI = sqrt(I);

        auto const numspecies = M.size();

        auto& ln_aw  = pzstate.ln_aw = 0.0;
        auto& OSMOT  = pzstate.phiw = 0.0;
        auto& LGAMMA = pzstate.ln_gamma = zeros(numspecies);

        real CSUM = 0.0;
        real BIGZ = (M * z.abs()).sum();
        real OSUM = M.sum() - M[iH2O];

        // If OSUM is zero, then solution is very dilluted - skip the rest and avoid division by zero when computing OSMOT
        if(OSUM == 0.0)
            return;

        // Compute the Debye-Huckel coefficient Aphi0 at (T, P)
        auto const Aphi0 = Aphi(T, P);

        // The b parameter of the Pitzer model
	    auto const B = 1.2;

        // The F term in the Pitzer model
	    real F = -Aphi0*(DI/(1 + B*DI) + 2.0*log(1.0 + B*DI)/B);

        // The osmotic coefficient of water in the Pitzer model
        OSMOT = -Aphi0*I*DI/(1 + B*DI);

        for(auto const& param : beta0)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            LGAMMA[i0] += M[i1] * 2.0 * param.value;
            LGAMMA[i1] += M[i0] * 2.0 * param.value;
            OSMOT += M[i0] * M[i1] * param.value;
        }

        for(auto const& [i, param] : enumerate(beta1))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            F += M[i0] * M[i1] * param.value * GP(alpha1[i] * DI)/I;
            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha1[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha1[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha1[i] * DI);
        }

        for(auto const& [i, param] : enumerate(beta2))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            F += M[i0] * M[i1] * param.value * GP(alpha2[i] * DI)/I;
            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha2[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha2[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha2[i] * DI);
        }

        for(auto const& param : Cphi)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            auto const aux = 2.0 * sqrt(abs(z[i0] * z[i1]));

            CSUM += M[i0] * M[i1] * param.value/aux;
            LGAMMA[i0] += M[i1] * BIGZ * param.value/aux;
            LGAMMA[i1] += M[i0] * BIGZ * param.value/aux;
            OSMOT += M[i0] * M[i1] * BIGZ * param.value/aux;
        }

        for(auto const& param : theta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            LGAMMA[i0] += 2.0 * M[i1] * param.value;
            LGAMMA[i1] += 2.0 * M[i0] * param.value;
            OSMOT += M[i0] * M[i1] * param.value;
        }

        for(auto const& [i0, i1] : thetaij)
        {
            auto const [etheta, ethetap] = computeThetaValuesInterpolation(I, DI, Aphi0, z[i0], z[i1]);

            F += M[i0] * M[i1] * ethetap;
            LGAMMA[i0] += 2.0 * M[i1] * etheta;
            LGAMMA[i1] += 2.0 * M[i0] * etheta;
            OSMOT += M[i0] * M[i1] * (etheta + I*ethetap);
        }

        for(auto const& param : psi)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        for(auto const& [i, param] : enumerate(lambda))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            auto const [clng0, clng1, cosm] = lambda_coeffs[i];

            LGAMMA[i0] += M[i1] * param.value * clng0;
            LGAMMA[i1] += M[i0] * param.value * clng1;
            OSMOT += M[i0] * M[i1] * param.value * cosm;
        }

        for(auto const& param : zeta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        for(auto const& [i, param] : enumerate(mu))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            auto const [clng0, clng1, clng2, cosm] = mu_coeffs[i];

            LGAMMA[i0] += M[i1] * M[i2] * param.value * clng0;
            LGAMMA[i1] += M[i0] * M[i2] * param.value * clng1;
            LGAMMA[i2] += M[i0] * M[i1] * param.value * clng2;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value * cosm;
        }

        for(auto const& param : eta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        // Finalise the calculation of the activity coefficient by adding the missing F and CSUM contributions
        for(auto i : icharged)
            LGAMMA[i] += z[i]*z[i]*F + abs(z[i])*CSUM;

        // Finalise the calculation of the activity coefficient by adding the missing F and CSUM contributions
        OSMOT = 1 + 2.0/OSUM * OSMOT;

        // Compute the activity of the water species
        ln_aw = -OSMOT * OSUM * Mw;
    }
};

} // namespace anonymous

auto createActivityModelExtendedUNIQUAC(SpeciesList const& species, ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModel
{
    // Create the aqueous solution
    AqueousMixture solution(species);

    // The index of water in the solution
    const Index iH2O = solution.indexWater();

    // Create the PitzerModel object responsible for evaluating the Pitzer model
    PitzerModel pzmodel(solution, params);

    // The PitzerState object that holds computed properties of the aqueous solution by the Pitzer model
    PitzerState pzstate;

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto aqstateptr = std::make_shared<AqueousMixtureState>();
    auto aqsolutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        // Evaluate the state of the aqueous solution
        auto const& aqstate = *aqstateptr = solution.state(T, P, x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // Export the aqueous solution and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = aqstateptr;
        props.extra["AqueousMixture"] = aqsolutionptr;

        // Evaluate the Pitzer activity model with given aqueous state
        pzmodel.evaluate(aqstate, pzstate);

        // The mole fraction of water
        auto const xw = x[iH2O];

        // Set the activity coefficients of the solutes
        props.ln_g = pzstate.ln_gamma;

        // Set the activities of the solutes
        props.ln_a = props.ln_g + log(aqstate.m);

        // Set the activitiy of water
        props.ln_a[iH2O] = pzstate.ln_aw;

        // Set the activity coefficient of water (mole fraction scale)
        props.ln_g[iH2O] = pzstate.ln_aw - log(xw);
    };

    return fn;
}

auto ActivityModelExtendedUNIQUAC() -> ActivityModelGenerator
{
    return ActivityModelExtendedUNIQUAC(Params::embedded("Pitzer.yaml"));
}

auto ActivityModelExtendedUNIQUAC(ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelExtendedUNIQUAC(species, params); };
}

auto ActivityModelExtendedUNIQUAC(Params const& params) -> ActivityModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ActivityModelParams"), "Expecting Pitzer activity model parameters in given Params object, but it lacks a `ActivityModelParams` section within which another section `Pitzer` should exist.");
    errorif(!data.at("ActivityModelParams").exists("Pitzer"), "Expecting Pitzer activity model parameters in given Params object, under the section `Pitzer`.");
    errorif(!data.at("ActivityModelParams").at("Pitzer").isDict(), "Expecting section `Pitzer` with Pitzer activity model parameters to be a dictionary.");

    ActivityModelParamsExtendedUNIQUAC pzparams =
        data["ActivityModelParams"]["Pitzer"].as<ActivityModelParamsExtendedUNIQUAC>();

    return ActivityModelExtendedUNIQUAC(pzparams);
}

} // namespace Reaktoro
