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

#include "ActivityModelPitzer.hpp"

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
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

using std::log;
using std::sqrt;
using std::exp;
using std::abs;
using std::round;

/// Auxiliary alias for ActivityModelParamsPitzer::CorrectionModel.
using PitzerParamCorrectionModel = ActivityModelParamsPitzer::CorrectionModel;

/// Return the Pitzer temperature-pressure correction model based on expression provided by PHREEQC v3.
auto createParamCorrectionModelPhreeqc(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;
    auto const Tr = 298.15;

    errorif(c.size() == 0, "Cannot create the Phreeqc temperature-dependent Pitzer parameter function with empty coefficients");

    if(c.size() == 1) return [=](real const& T, real const& Pbar) { return c[0]; };
    if(c.size() == 2) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr); };
    if(c.size() == 3) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr); };
    if(c.size() == 4) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr); };
    if(c.size() == 5) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr); };
    if(c.size() == 6) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr) + c[5]*(1.0/(T*T) - 1.0/(Tr*Tr)); };

    errorif(true, "Cannot create the Phreeqc temperature-dependent Pitzer parameter function with given coefficients: ", str(c));
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by He and Morse (1993) (doi: 10.1016/0016-7037(93)90137-L).
auto createParamCorrectionModelHeMorse1993(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `HeMorse1993` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Dai et al. (2013) (doi: 10.2118/164045-ms).
auto createParamCorrectionModelDai2013(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Dai2013` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Dai et al. (2014) (doi: 10.2118/169786-ms).
auto createParamCorrectionModelDai2014(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Dai2014` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Christov and Møller (2004) (see Table 5 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelChristovMoller2004(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `ChristovMoller2004` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Holmes et al. (1987) (see Table 6 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelHolmes1987(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Holmes1987` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Pitzer et al. (1984) (see Table 7 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPitzer1984(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Pitzer1984` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Palaban and Pitzer (1987) (see Table 8 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPalabanPitzer1987(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `PalabanPitzer1987` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Polya et al. (2001) (see Table 10 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPolya2001(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Polya2001` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Lia and Duan (2007) (see Table 12 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelLiDuan2007(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `LiDuan2007` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model with given coefficients and model option.
auto createParamCorrectionModel(Vec<Param> const& coefficients, PitzerParamCorrectionModel const& option) -> Fn<real(real const&, real const&)>
{
    switch(option)
    {
        case PitzerParamCorrectionModel::Phreeqc:            return createParamCorrectionModelPhreeqc(coefficients);
        case PitzerParamCorrectionModel::HeMorse1993:        return createParamCorrectionModelHeMorse1993(coefficients);
        case PitzerParamCorrectionModel::Dai2013:            return createParamCorrectionModelDai2013(coefficients);
        case PitzerParamCorrectionModel::Dai2014:            return createParamCorrectionModelDai2014(coefficients);
        case PitzerParamCorrectionModel::ChristovMoller2004: return createParamCorrectionModelChristovMoller2004(coefficients);
        case PitzerParamCorrectionModel::Holmes1987:         return createParamCorrectionModelHolmes1987(coefficients);
        case PitzerParamCorrectionModel::Pitzer1984:         return createParamCorrectionModelPitzer1984(coefficients);
        case PitzerParamCorrectionModel::PalabanPitzer1987:  return createParamCorrectionModelPalabanPitzer1987(coefficients);
        case PitzerParamCorrectionModel::Polya2001:          return createParamCorrectionModelPolya2001(coefficients);
        case PitzerParamCorrectionModel::LiDuan2007:         return createParamCorrectionModelLiDuan2007(coefficients);
        default:                                             return createParamCorrectionModelPhreeqc(coefficients);
    }
}

/// Used to represent an interaction parameter in the Pitzer activity model.
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

/// Auxiliary alias for ActivityModelParamsPitzer::InteractionParamAttribs.
using PitzerInteractionParamAttribs = ActivityModelParamsPitzer::InteractionParamAttribs;

/// Return a given list of species indices sorted in ascending order of electric charge.
auto sortedSpeciesIndicesByCharge(SpeciesList const& specieslist, Indices const& ispecies) -> Indices
{
    // Lambda function that returns the charge of the ith species in `specieslist`
    auto chargefn = [&](auto i) { return specieslist[i].charge(); };

    // Lambda function that returns true if two species indices are in ascending order of electric charge.
    auto lesschargefn = [&](Index ispecies1, Index ispecies2) { return charge(ispecies1) < charge(ispecies2); };

    return sortedfn(ispecies, lesschargefn);
}

/// Convert a PitzerInteractionParamAttribs object to a PitzerParam one for a binary interaction parameter.
auto createPitzerParamBinary(SpeciesList const& specieslist, PitzerInteractionParamAttribs const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 2, "Expecting two chemical formulas when processing a Pitzer binary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);

    if(ispecies1 >= specieslist.size() || ispecies2 >= specieslist.size())
        return {}; // specie1 and/or species2 are not in the list of species; return empty PitzerParam object.

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.coefficients, attribs.model) };
}

/// Convert a PitzerInteractionParamAttribs object to a PitzerParam one for a ternary interaction parameter.
auto createPitzerParamTernary(SpeciesList const& specieslist, PitzerInteractionParamAttribs const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 3, "Expecting two chemical formulas when processing a Pitzer ternary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);
    auto const ispecies3 = specieslist.findWithFormula(attribs.formulas[2]);

    if(ispecies1 >= specieslist.size() || ispecies2 >= specieslist.size() || ispecies3 >= specieslist.size())
        return {}; // specie1 and/or species2 and/or species3 are not in the list of species; return empty PitzerParam object.

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2, ispecies3});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.coefficients, attribs.model) };
}

/// Auxiliary alias for ActivityModelParamsPitzer::AlphaParamAttribs
using PitzerAlphaParamAttribs = ActivityModelParamsPitzer::AlphaParamAttribs;

/// Return the default value for \eq{alpha_1} parameter according to that used in PHREEQC v3 (see file pitzer.cpp under comment "Set alpha values").
auto determineDefaultAlpha1(ChemicalFormula const& formula1, ChemicalFormula const& formula2) -> Param
{
    const auto z1 = round(abs(formula1.charge()));
    const auto z2 = round(abs(formula2.charge()));
    return (z1 == 2 && z2 == 2) ? 1.4 : 2.0; // alpha1 = 2.0 for all electrolytes except the 2-2 electrolytes, for which alpha1 = 1.4 (Pitzer and Silvester, 1978; Plummer et al., 1988, page 4)
};

/// Return the default value for \eq{alpha_2} parameter according to that used in PHREEQC v3 (see file pitzer.cpp under comment "Set alpha values").
auto determineDefaultAlpha2(ChemicalFormula const& formula1, ChemicalFormula const& formula2) -> Param
{
    const auto z1 = round(abs(formula1.charge()));
    const auto z2 = round(abs(formula2.charge()));
    return (z1 == 1 && z2 == 1) || (z1 == 2 && z2 == 2) ? 12.0 : 50.0; // alpha2 = 12.0 for 1-1 and 2-2 electrolytes, alpha2 = 50.0 for 3-2 and 4-2 electrolytes (Pitzer and Silvester, 1978; Plummer et al., 1988, page 4)
};

/// Return the index of the entry in `alphas` containing a pair of formulas equivalent to the pair `formula1` and `formula2`.
auto findAlphaEntry(ChemicalFormula const& formula1, ChemicalFormula const& formula2, Vec<PitzerAlphaParamAttribs> const& alphas) -> Index
{
    return indexfn(alphas, RKT_LAMBDA(x,
        ( formula1.equivalent(x.formula1) && formula2.equivalent(x.formula2) ) ||
        ( formula1.equivalent(x.formula2) && formula2.equivalent(x.formula1) )
    ));
}

/// Return the value for \eq{alpha_1} parameter (either user-specified or determined from default value according to that used in PHREEQC v3).
auto determineAlpha1(ChemicalFormula const& formula1, ChemicalFormula const& formula2, Vec<PitzerAlphaParamAttribs> const& alphas) -> Param
{
    if(const auto idx = findAlphaEntry(formula1, formula2, alphas); idx < alphas.size())
        return alphas[idx].value;
    return determineDefaultAlpha1(formula1, formula2);
}

/// Return the value for \eq{alpha_2} parameter (either user-specified or determined from default value according to that used in PHREEQC v3).
auto determineAlpha2(ChemicalFormula const& formula1, ChemicalFormula const& formula2, Vec<PitzerAlphaParamAttribs> const& alphas) -> Param
{
    if(const auto idx = findAlphaEntry(formula1, formula2, alphas); idx < alphas.size())
        return alphas[idx].value;
    return determineDefaultAlpha2(formula1, formula2);
}

/// Determine the coefficients multiplying the terms where the lambda Pitzer parameter is involved.
/// Check equations (34-37) of Clegg (1991) (Activity Coefficients in Natural Waters).
auto determineLambdaCoeffs(double z1, double z2, Index i1, Index i2) -> Tuple<double, double, double>
{
    assert(z1 <= z2 && "Expecting species associated to binary lambda Pitzer parameter to be ordered in ascending order of charge.");

    const bool NxNx = z1 == 0.0 && z2 == 0.0 && i1 == i2;
    if(NxNx)
        return {1.0, 1.0, 0.5}; // (37.i => 2/2), (37.i => 2/2), (34.iv => 1/2)

    const bool NxNy = z1 == 0.0 && z2 == 0.0 && i1 != i2;
    if(NxNy)
        return {2.0, 2.0, 1.0}; // (37.i => 2), (37.i => 2), (34.v => 1)

    const bool NxCx = z1 == 0.0 && z2  > 0.0;
    if(NxCx)
        return {2.0, 2.0, 1.0}; // (37.iv => 2), (35.iv => 2), (34.vii => 1)

    const bool AxNx = z1  < 0.0 && z2 == 0.0;
    if(AxNx)
        return {2.0, 2.0, 1.0}; // (36.iv => 2), (37.iv => 2), (34.viii => 1)

    errorif(true, "Could not determine mu factor in osmotic coefficient equation with given triad of species with charges (", z1, ", ", z2, ").");

    return {};
}

/// Determine the coefficients multiplying the terms where the mu Pitzer parameter is involved.
/// Check equations (34-37) of Clegg (1991) (Activity Coefficients in Natural Waters).
auto determineMuCoeffs(double z1, double z2, double z3, Index i1, Index i2, Index i3) -> Tuple<double, double, double, double>
{
    assert(z1 <= z2 && z2 <= z3 && "Expecting species associated to ternary mu Pitzer parameter to be ordered in ascending order of charge.");

    const bool NNN = z1 == 0.0 && z2 == 0.0 && z3 == 0.0;
    const bool NNC = z1 == 0.0 && z2 == 0.0 && z3  > 0.0;
    const bool ANN = z1  < 0.0 && z2 == 0.0 && z3 == 0.0;

    const bool NxNxNx = NNN && i1 == i2 && i2 == i3;
    if(NxNxNx)
        return {1.0, 1.0, 1.0, 1.0}; // (37.i => 3/3), (37.i => 3/3), (37.i => 3/3), (34.iv => 1)

    const bool NxNxNy = NNN && i1 == i2 && i2 != i3;
    if(NxNxNy)
        return {3.0, 3.0, 3.0, 3.0}; // (37.ii => 6/2), (37.ii => 6/2), (37.i => 3), (34.v => 3)

    const bool NyNxNx = NNN && i1 != i2 && i2 == i3;
    if(NyNxNx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.i => 3), (37.ii => 6/2), (37.ii => 6/2), (34.v => 3)

    const bool NxNyNx = NNN && i1 != i2 && i1 == i3;
    if(NxNyNx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.ii => 6/2), (37.i => 3), (37.ii => 6/2), (34.v => 3)

    const bool NxNyNz = NNN && i1 != i2 && i1 != i3;
    if(NxNyNz)
        return {6.0, 6.0, 6.0, 6.0}; // (37.iii => 6), (37.iii => 6), (37.iii => 6), (34.vi => 6)

    const bool NxNxCx = NNC && i1 == i2;
    if(NxNxCx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.vii => 6/2), (37.vii => 6/2), (35.iv => 3), (34.vii => 3)

    const bool AxNxNx = ANN && i2 == i3;
    if(AxNxNx)
        return {3.0, 3.0, 3.0, 3.0}; // (36.iv => 3), (37.vii => 6/2), (37.vii => 6/2), (34.viii => 3)

    const bool NxNyCx = NNC && i1 != i2;
    if(NxNyCx)
        return {6.0, 6.0, 6.0, 6.0}; // (37.vii => 6), (37.vii => 6), (35.v => 6), (34.xii => 6)

    const bool AxNxNy = ANN && i2 != i3;
    if(AxNxNy)
        return {6.0, 6.0, 6.0, 6.0}; // (36.v => 6), (37.vii => 6/2), (37.vii => 6/2), (34.xiii => 6)

    errorif(true, "Could not determine mu factor in osmotic coefficient equation with given triad of species with charges (", z1, ", ", z2, ", ", z3, ").");

    return {};
}

auto G(real const& x) -> real
{
    return x == 0.0 ? x : 2.0*(1.0 - (1.0 + x)*exp(-x))/(x*x);
}

auto GP(real const& x) -> real
{
    return x == 0.0 ? x : -2.0*(1.0 - (1.0 + x + 0.5*x*x)*exp(-x))/(x*x);
}

auto computeThetaValues(real const& I, real const& sqrtI, real const& Aphi, double zi, double zj) -> Pair<real, real>
{
    const auto aux = 6.0*Aphi*sqrtI;

    const auto xij = zi*zj*aux;
    const auto xii = zi*zi*aux;
    const auto xjj = zj*zj*aux;

    const auto J0ij = J0(xij);
    const auto J0ii = J0(xii);
    const auto J0jj = J0(xjj);

    const auto J1ij = J1(xij);
    const auto J1ii = J1(xii);
    const auto J1jj = J1(xjj);

    const auto thetaE = zi*zj/(4*I) * (J0ij - 0.5*J0ii - 0.5*J0jj);
    const auto thetaEP = zi*zj/(8*I*I) * (J1ij - 0.5*J1ii - 0.5*J1jj) - thetaE/I;

    return { thetaE, thetaEP };
}

/// The auxiliary type used to store computed values from the Pitzer activity model implemented by @ref Pitzer.
struct PitzerState
{
    ArrayXr phiw; ///< The osmotic coefficient of water.
    ArrayXr ln_g; ///< The activity coefficients of the aqueous species (natural log).
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

    using Tuples3d = Tuples<double, double, double>; ///< Auxiliary type for a tuple of 3 double values.
    using Tuples4d = Tuples<double, double, double, double>; ///< Auxiliary type for a tuple of 4 double values.

    Tuples3d lambda_coeffs; ///< The coefficients multiplying the terms where the lambda Pitzer parameter is involved.
    Tuples4d mu_coeffs;     ///< The coefficients multiplying the terms where the mu Pitzer parameter is involved.

    ArrayXr thetaE;  ///< The current values of the parameters \eq{^{E}\theta_{ij}(I)} associated to the \eq{\theta_{ij}} parameters.
    ArrayXr thetaEP; ///< The current values of the parameters \eq{^{E}\theta_{ij}^{\prime}(I)} associated to the \eq{\theta_{ij}} parameters.

    BilinearInterpolator Aphi; ///< The parameter \eq{A^\phi(T, P)} in the Pitzer model.

    /// Construct a default Pitzer object.
    PitzerModel()
    {}

    /// Construct a Pitzer object with given list of species in the aqueous solution and the parameters for the Pitzer activity model.
    // PitzerModel(SpeciesList const& specieslist, ActivityModelParamsPitzer const& params)
    PitzerModel(AqueousMixture const& solution, ActivityModelParamsPitzer const& params)
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

        auto const& z = solution.charges();

        for(auto const& entry : lambda)
        {
            const auto i1 = entry.ispecies[0];
            const auto i2 = entry.ispecies[1];
            lambda_coeffs.push_back(determineLambdaCoeffs(z[i1], z[i2], i1, i2));
        }

        for(auto const& entry : mu)
        {
            const auto i1 = entry.ispecies[0];
            const auto i2 = entry.ispecies[1];
            const auto i3 = entry.ispecies[2];
            mu_coeffs.push_back(determineMuCoeffs(z[i1], z[i2], z[i3], i1, i2, i3));
        }

        Aphi = BilinearInterpolator(Aphi_temperatures, Aphi_pressures, Aphi_data);
    }

    auto evaluate(AqueousMixtureState const& aqstate, PitzerState& pzstate)
    {
        auto const& T = aqstate.T;
        auto const& P = aqstate.P;
        auto const& M = aqstate.m;
        auto const& I = aqstate.Ie;
        auto const& z = solution.charges();

        auto const& iH2O = solution.indexWater();

        auto const DI = sqrt(I);

        auto const numspecies = M.size();

        auto& OSMOT  = pzstate.phiw = 0.0;
        auto& LGAMMA = pzstate.ln_g = zeros(numspecies);

        real CSUM = 0.0;

        real BIGZ = (M * z.abs()).sum();
        real OSUM = M.sum() - M[iH2O];

        real Aphi0 = Aphi(T, P);

	    const auto B = 1.2;

	    real F = -Aphi0*(DI/(1 + B*DI) + 2.0*log(1.0 + B*DI)/B);

        for(auto const& param : beta0)
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];

			LGAMMA[i0] += M[i1] * 2.0 * param.value;
			LGAMMA[i1] += M[i0] * 2.0 * param.value;
			OSMOT += M[i0] * M[i1] * param.value;
        }

        for(auto const& [i, param] : enumerate(beta1))
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];

            F += M[i0] * M[i1] * param.value * GP(alpha1[i] * DI) / I;
            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha1[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha1[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha1[i] * DI);
        }

        for(auto const& [i, param] : enumerate(beta2))
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];

            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha2[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha2[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha2[i] * DI);
        }

        for(auto const& param : Cphi)
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];
            const auto z0 = z[i0];
            const auto z1 = z[i1];
            const auto aux = 2.0 * sqrt(abs(z0 * z1));

            CSUM += M[i0] * M[i1] * param.value / aux;
            LGAMMA[i0] += M[i1] * BIGZ * param.value / aux;
            LGAMMA[i1] += M[i0] * BIGZ * param.value / aux;
            OSMOT += M[i0] * M[i1] * BIGZ * param.value / aux;
        }

        for(auto const& param : theta)
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];
            const auto z0 = z[i0];
            const auto z1 = z[i1];

            LGAMMA[i0] += 2.0 * M[i1] * param.value;
            LGAMMA[i1] += 2.0 * M[i0] * param.value;
            OSMOT += M[i0] * M[i1] * param.value;

            const auto [etheta, ethetap] = computeThetaValues(I, DI, Aphi0, z0, z1);

            F += M[i0] * M[i1] * ethetap;
            LGAMMA[i0] += 2.0 * M[i1] * etheta;
            LGAMMA[i1] += 2.0 * M[i0] * etheta;
            OSMOT += M[i0] * M[i1] * (etheta + I * ethetap);
        }

        for(auto const& param : psi)
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];
            const auto i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        for(auto const& param : lambda)
        {
            const auto i0 = param.ispecies[0];
            const auto i1 = param.ispecies[1];

            LGAMMA[i0] += M[i1] * param.value * pitz_params[i]->ln_coef[0];
            LGAMMA[i1] += M[i0] * param.value * pitz_params[i]->ln_coef[1];
            OSMOT += M[i0] * M[i1] * param.value * pitz_params[i]->os_coef;
        }
    }
};


template<typename T> using Table1D = Vec<T>;
template<typename T> using Table2D = Vec<Vec<T>>;
template<typename T> using Table3D = Vec<Vec<Vec<T>>>;

template<typename T>
auto initTable1D(unsigned dim1) -> Table1D<T>
{
    return Table1D<T>(dim1);
}

template<typename T>
auto initTable2D(unsigned dim1, unsigned dim2) -> Table2D<T>
{
    return Table2D<T>(dim1, initTable1D<T>(dim2));
}

template<typename T>
auto initTable3D(unsigned dim1, unsigned dim2, unsigned dim3) -> Table3D<T>
{
    return Table3D<T>(dim1, initTable2D<T>(dim2, dim3));
}

const Pairs<Vec<ChemicalFormula>, Vec<real>> beta0_data =
{
    { {"Ba++"     , "Br-"        }, {  0.31455, -3.3825e-4                                  } },
    { {"Ba++"     , "Cl-"        }, {  0.26280,  6.4050e-4                                  } },
    { {"Ba++"     , "OH-"        }, {  0.17175                                              } },
    { {"Ca++"     , "Br-"        }, {  0.38160, -5.2275e-4                                  } },
    { {"Ca++"     , "Cl-"        }, {  0.31590, -1.7250e-4                                  } },
    { {"Ca++"     , "HCO3-"      }, {  0.40000                                              } },
    { {"Ca++"     , "HSO4-"      }, {  0.21450                                              } },
    { {"Ca++"     , "OH-"        }, { -0.17470                                              } },
    { {"Ca++"     , "SO4--"      }, {  0.20000                                              } },
    { {"CaB(OH)4+", "Cl-"        }, {  0.12000                                              } },
    { {"Fe++"     , "Cl-"        }, {  0.33593                                              } },
    { {"Fe++"     , "HSO4-"      }, {  0.42730                                              } },
    { {"Fe++"     , "SO4--"      }, {  0.25680                                              } },
    { {"H+"       , "Br-"        }, {  0.19600, -2.0490e-4                                  } },
    { {"H+"       , "Cl-"        }, {  0.17750, -3.0810e-4                                  } },
    { {"H+"       , "HSO4-"      }, {  0.20650                                              } },
    { {"H+"       , "SO4--"      }, {  0.02980                                              } },
    { {"K+"       , "B(OH)4-"    }, {  0.03500                                              } },
    { {"K+"       , "B3O3(OH)4-" }, { -0.13000                                              } },
    { {"K+"       , "B4O5(OH)4--"}, { -0.02200                                              } },
    { {"K+"       , "Br-"        }, {  0.05690,  7.3900e-4                                  } },
    { {"K+"       , "Cl-"        }, {  0.04835,  5.7940e-4                                  } },
    { {"K+"       , "CO3--"      }, {  0.14880,  1.7880e-3                                  } },
    { {"K+"       , "HCO3-"      }, {  0.02960,  9.9600e-4                                  } },
    { {"K+"       , "HSO4-"      }, { -0.00030                                              } },
    { {"K+"       , "OH-"        }, {  0.12980                                              } },
    { {"K+"       , "SO4--"      }, {  0.04995,  1.4400e-3                                  } },
    { {"Li+"      , "Br-"        }, {  0.17480, -1.8190e-4                                  } },
    { {"Li+"      , "Cl-"        }, {  0.14940, -1.6850e-4                                  } },
    { {"Li+"      , "OH-"        }, {  0.01500                                              } },
    { {"Li+"      , "SO4--"      }, {  0.13628,  5.0550e-4                                  } },
    { {"Mg++"     , "Br-"        }, {  0.43270, -5.6250e-5                                  } },
    { {"Mg++"     , "Cl-"        }, {  0.35235, -1.9430e-4                                  } },
    { {"Mg++"     , "HCO3-"      }, {  0.32900                                              } },
    { {"Mg++"     , "HSO4-"      }, {  0.47460                                              } },
    { {"Mg++"     , "SO4--"      }, {  0.22100, -6.9000e-4                                  } },
    { {"MgB(OH)4+", "Cl-"        }, {  0.16000                                              } },
    { {"MgOH+"    , "Cl-"        }, { -0.10000                                              } },
    { {"Mn++"     , "Cl-"        }, {  0.32723                                              } },
    { {"Mn++"     , "SO4--"      }, {  0.20650                                              } },
    { {"Na+"      , "B(OH)4-"    }, { -0.04270                                              } },
    { {"Na+"      , "B3O3(OH)4-" }, { -0.05600                                              } },
    { {"Na+"      , "B4O5(OH)4--"}, { -0.11000                                              } },
    { {"Na+"      , "Br-"        }, {  0.09730,  7.6920e-4                                  } },
    { {"Na+"      , "Cl-"        }, {  0.07650, -777.03    , -4.4706, 8.9460e-3, -3.3158e-6 } },
    { {"Na+"      , "CO3--"      }, {  0.03990,  1.7900e-3                                  } },
    { {"Na+"      , "HCO3-"      }, {  0.02770,  1.0000e-3                                  } },
    { {"Na+"      , "HSO4-"      }, {  0.04540                                              } },
    { {"Na+"      , "OH-"        }, {  0.08640,  7.0000e-4                                  } },
    { {"Na+"      , "SO4--"      }, {  0.01958,  2.3670e-3                                  } },
    { {"Sr++"     , "Br-"        }, {  0.33113, -3.2775e-4                                  } },
    { {"Sr++"     , "Cl-"        }, {  0.28580,  7.1700e-4                                  } },
    { {"Sr++"     , "HCO3-"      }, {  0.12000                                              } },
    { {"Sr++"     , "SO4--"      }, {  0.20000, -2.9000e-3                                  } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> beta1_data =
{
    { {"Ba++" , "Br-"        }, {  1.56975, 6.7800e-3                                  } },
    { {"Ba++" , "Cl-"        }, {  1.49625, 3.2325e-3                                  } },
    { {"Ba++" , "OH-"        }, {  1.20000                                             } },
    { {"Ca++" , "Br-"        }, {  1.61300, 6.0375e-3                                  } },
    { {"Ca++" , "Cl-"        }, {  1.61400, 3.9000e-3                                  } },
    { {"Ca++" , "HCO3-"      }, {  2.97700                                             } },
    { {"Ca++" , "HSO4-"      }, {  2.53000                                             } },
    { {"Ca++" , "OH-"        }, { -0.23030                                             } },
    { {"Ca++" , "SO4--"      }, {  3.19730, 5.4600e-2                                  } },
    { {"Fe++" , "Cl-"        }, {  1.53225                                             } },
    { {"Fe++" , "HSO4-"      }, {  3.48000                                             } },
    { {"Fe++" , "SO4--"      }, {  3.06300                                             } },
    { {"H+"   , "Br-"        }, {  0.35640, 4.4670e-4                                  } },
    { {"H+"   , "Cl-"        }, {  0.29450, 1.4190e-4                                  } },
    { {"H+"   , "HSO4-"      }, {  0.55560                                             } },
    { {"K+"   , "B(OH)4-"    }, {  0.14000                                             } },
    { {"K+"   , "Br-"        }, {  0.22120, 1.7400e-3                                  } },
    { {"K+"   , "Cl-"        }, {  0.21220, 1.0710e-3                                  } },
    { {"K+"   , "CO3--"      }, {  1.43000, 2.0510e-3                                  } },
    { {"K+"   , "HCO3-"      }, { -0.01300, 1.1040e-3                                  } },
    { {"K+"   , "HSO4-"      }, {  0.17350                                             } },
    { {"K+"   , "OH-"        }, {  0.32000                                             } },
    { {"K+"   , "SO4--"      }, {  0.77930, 6.6975e-3                                  } },
    { {"Li+"  , "Br-"        }, {  0.25470, 6.6360e-4                                  } },
    { {"Li+"  , "Cl-"        }, {  0.30740, 5.3660e-4                                  } },
    { {"Li+"  , "OH-"        }, {  0.14000                                             } },
    { {"Li+"  , "SO4--"      }, {  1.27050, 1.4100e-3                                  } },
    { {"Mg++" , "Br-"        }, {  1.75300, 3.8625e-3                                  } },
    { {"Mg++" , "Cl-"        }, {  1.68150, 3.6525e-3                                  } },
    { {"Mg++" , "HCO3-"      }, {  0.60720                                             } },
    { {"Mg++" , "HSO4-"      }, {  1.72900                                             } },
    { {"Mg++" , "SO4--"      }, {  3.34300, 1.5300e-2                                  } },
    { {"MgOH+", "Cl-"        }, {  1.65800                                             } },
    { {"Mn++" , "Cl-"        }, {  1.55025                                             } },
    { {"Mn++" , "SO4--"      }, {  2.95110                                             } },
    { {"Na+"  , "B(OH)4-"    }, {  0.08900                                             } },
    { {"Na+"  , "B3O3(OH)4-" }, { -0.91000                                             } },
    { {"Na+"  , "B4O5(OH)4--"}, { -0.40000                                             } },
    { {"Na+"  , "Br-"        }, {  0.27910, 1.0790e-3                                  } },
    { {"Na+"  , "Cl-"        }, {  0.26640, 0.0000e+0, 0.0000e+0, 6.1608e-5, 1.0715e-6 } },
    { {"Na+"  , "CO3--"      }, {  1.38900, 2.0500e-3                                  } },
    { {"Na+"  , "HCO3-"      }, {  0.04110, 1.1000e-3                                  } },
    { {"Na+"  , "HSO4-"      }, {  0.39800                                             } },
    { {"Na+"  , "OH-"        }, {  0.25300, 1.3400e-4                                  } },
    { {"Na+"  , "SO4--"      }, {  1.11300, 5.6325e-3                                  } },
    { {"Sr++" , "Br-"        }, {  1.71150, 6.5325e-3                                  } },
    { {"Sr++" , "Cl-"        }, {  1.66700, 2.8425e-3                                  } },
    { {"Sr++" , "SO4--"      }, {  3.19730, 2.7000e-2                                  } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> beta2_data =
{
    { {"Ca++",  "OH-"   }, { -5.720         }},
    { {"Ca++",  "SO4--" }, { -54.24, -0.516 }},
    { {"Fe++",  "SO4--" }, { -42.00         }},
    { {"Mg++",  "SO4--" }, { -37.23, -0.253 }},
    { {"Mn++",  "SO4--" }, { -40.00         }},
    { {"Sr++",  "SO4--" }, { -54.24, -0.420 }}
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> Cphi_data =
{
    { {"Ba++", "Br-"    }, { -0.01595760                                            } },
    { {"Ba++", "Cl-"    }, { -0.01937820, -1.53796e-4                               } },
    { {"Ca++", "Br-"    }, { -0.00257000                                            } },
    { {"Ca++", "Cl-"    }, { -0.00034000                                            } },
    { {"Fe++", "Cl-"    }, { -0.00860725                                            } },
    { {"Fe++", "SO4--"  }, {  0.02090000                                            } },
    { {"H+"  , "Br-"    }, {  0.00827000, -5.68500e-5                               } },
    { {"H+"  , "Cl-"    }, {  0.00080000,  6.21300e-5                               } },
    { {"H+"  , "SO4--"  }, {  0.04380000                                            } },
    { {"K+"  , "Br-"    }, { -0.00180000, -7.00400e-5                               } },
    { {"K+"  , "Cl-"    }, { -0.00084000, -5.09500e-5                               } },
    { {"K+"  , "CO3--"  }, { -0.00150000                                            } },
    { {"K+"  , "HCO3-"  }, { -0.00800000                                            } },
    { {"K+"  , "OH-"    }, {  0.00410000                                            } },
    { {"Li+" , "Br-"    }, {  0.00530000, -2.81300e-5                               } },
    { {"Li+" , "Cl-"    }, {  0.00359000, -4.52000e-5                               } },
    { {"Li+" , "SO4--"  }, { -0.00399338, -2.33345e-4                               } },
    { {"Mg++", "Br-"    }, {  0.00312000                                            } },
    { {"Mg++", "Cl-"    }, {  0.00519000, -1.64933e-4                               } },
    { {"Mg++", "SO4--"  }, {  0.02500000,  5.23000e-4                               } },
    { {"Mn++", "Cl-"    }, { -0.02049720                                            } },
    { {"Mn++", "SO4--"  }, {  0.01636000                                            } },
    { {"Na+" , "B(OH)4-"}, {  0.01140000                                            } },
    { {"Na+" , "Br-"    }, {  0.00116000, -9.30000e-5                               } },
    { {"Na+" , "Cl-"    }, {  0.00127000,  33.317, 0.09421, -4.65500e-5, 0.00000e+0 } },
    { {"Na+" , "CO3--"  }, {  0.00440000                                            } },
    { {"Na+" , "OH-"    }, {  0.00440000, -1.89400e-4                               } },
    { {"Na+" , "SO4--"  }, {  0.00497000, -4.87904e-4                               } },
    { {"Sr++", "Br-"    }, {  0.00122506                                            } },
    { {"Sr++", "Cl-"    }, { -0.00130000                                            } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> theta_data =
{
    { {"B(OH)4-"    , "Cl-"  }, { -0.065 } },
    { {"B(OH)4-"    , "SO4--"}, { -0.012 } },
    { {"B3O3(OH)4-" , "Cl-"  }, {  0.120 } },
    { {"B3O3(OH)4-" , "HCO3-"}, { -0.100 } },
    { {"B3O3(OH)4-" , "SO4--"}, {  0.100 } },
    { {"B4O5(OH)4--", "Cl-"  }, {  0.074 } },
    { {"B4O5(OH)4--", "HCO3-"}, { -0.087 } },
    { {"B4O5(OH)4--", "SO4--"}, {  0.120 } },
    { {"Br-"        , "OH-"  }, { -0.065 } },
    { {"Ca++"       , "H+"   }, {  0.092 } },
    { {"Ca++"       , "K+"   }, {  0.032 } },
    { {"Ca++"       , "Mg++" }, {  0.007 } },
    { {"Ca++"       , "Na+"  }, {  0.070 } },
    { {"Cl-"        , "CO3--"}, { -0.020 } },
    { {"Cl-"        , "HCO3-"}, {  0.030 } },
    { {"Cl-"        , "HSO4-"}, { -0.006 } },
    { {"Cl-"        , "OH-"  }, { -0.050 } },
    { {"Cl-"        , "SO4--"}, {  0.020 } },
    { {"CO3--"      , "OH-"  }, {  0.100 } },
    { {"CO3--"      , "SO4--"}, {  0.020 } },
    { {"H+"         , "K+"   }, {  0.005 } },
    { {"H+"         , "Mg++" }, {  0.100 } },
    { {"H+"         , "Na+"  }, {  0.036 } },
    { {"HCO3-"      , "CO3--"}, { -0.040 } },
    { {"HCO3-"      , "SO4--"}, {  0.010 } },
    { {"K+"         , "Na+"  }, { -0.012 } },
    { {"Mg++"       , "Na+"  }, {  0.070 } },
    { {"Na+"        , "Sr++" }, {  0.051 } },
    { {"OH-"        , "SO4--"}, { -0.013 } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> psi_data =
{
    { {"B(OH)4-",    "Cl-"  , "Na+"  }, { -0.0073 } },
    { {"B3O3(OH)4-", "Cl-"  , "Na+"  }, { -0.0240 } },
    { {"B4O5(OH)4--","Cl-"  , "Na+"  }, {  0.0260 } },
    { {"Br-",        "K+"   , "Na+"  }, { -0.0022 } },
    { {"Br-",        "K+"   , "OH-"  }, { -0.0140 } },
    { {"Br-",        "Na+"  , "H+"   }, { -0.0120 } },
    { {"Br-",        "Na+"  , "OH-"  }, { -0.0180 } },
    { {"Ca++",       "Cl-"  , "H+"   }, { -0.0150 } },
    { {"Ca++",       "Cl-"  , "K+"   }, { -0.0250 } },
    { {"Ca++",       "Cl-"  , "Mg++" }, { -0.0120 } },
    { {"Ca++",       "Cl-"  , "Na+"  }, { -0.0070 } },
    { {"Ca++",       "Cl-"  , "OH-"  }, { -0.0250 } },
    { {"Ca++",       "Cl-"  , "SO4--"}, { -0.0180 } },
    { {"Ca++",       "Mg++" , "SO4--"}, {  0.0240 } },
    { {"Ca++",       "Na+"  , "SO4--"}, { -0.0550 } },
    { {"Cl-",        "Br-"  , "K+"   }, {  0.0000 } },
    { {"Cl-",        "CO3--", "K+"   }, {  0.0040 } },
    { {"Cl-",        "CO3--", "Na+"  }, {  0.0085 } },
    { {"Cl-",        "H+"   , "K+"   }, { -0.0110 } },
    { {"Cl-",        "H+"   , "Mg++" }, { -0.0110 } },
    { {"Cl-",        "H+"   , "Na+"  }, { -0.0040 } },
    { {"Cl-",        "HCO3-", "Mg++" }, { -0.0960 } },
    { {"Cl-",        "HCO3-", "Na+"  }, { -0.0150 } },
    { {"Cl-",        "HSO4-", "H+"   }, {  0.0130 } },
    { {"Cl-",        "HSO4-", "Na+"  }, { -0.0060 } },
    { {"Cl-",        "K+"   , "Mg++" }, { -0.0220 } },
    { {"Cl-",        "K+"   , "Na+"  }, { -0.0018 } },
    { {"Cl-",        "K+"   , "OH-"  }, { -0.0060 } },
    { {"Cl-",        "Mg++" , "MgOH+"}, {  0.0280 } },
    { {"Cl-",        "Mg++" , "Na+"  }, { -0.0120 } },
    { {"Cl-",        "Mg++" , "SO4--"}, { -0.0040 } },
    { {"Cl-",        "Na+"  , "OH-"  }, { -0.0060 } },
    { {"Cl-",        "Na+"  , "SO4--"}, {  0.0014 } },
    { {"Cl-",        "Na+"  , "Sr++" }, { -0.0021 } },
    { {"CO3--",      "HCO3-", "K+"   }, {  0.0120 } },
    { {"CO3--",      "HCO3-", "Na+"  }, {  0.0020 } },
    { {"CO3--",      "K+"   , "Na+"  }, {  0.0030 } },
    { {"CO3--",      "K+"   , "OH-"  }, { -0.0100 } },
    { {"CO3--",      "K+"   , "SO4--"}, { -0.0090 } },
    { {"CO3--",      "Na+"  , "OH-"  }, { -0.0170 } },
    { {"CO3--",      "Na+"  , "SO4--"}, { -0.0050 } },
    { {"H+",         "HSO4-", "K+"   }, { -0.0265 } },
    { {"H+",         "HSO4-", "Mg++" }, { -0.0178 } },
    { {"H+",         "HSO4-", "Na+"  }, { -0.0129 } },
    { {"H+",         "K+"   , "Br-"  }, { -0.0210 } },
    { {"H+",         "K+"   , "SO4--"}, {  0.1970 } },
    { {"HCO3-",      "K+"   , "Na+"  }, { -0.0030 } },
    { {"HCO3-",      "Mg++" , "SO4--"}, { -0.1610 } },
    { {"HCO3-",      "Na+"  , "SO4--"}, { -0.0050 } },
    { {"HSO4-",      "K+"   , "SO4--"}, { -0.0677 } },
    { {"HSO4-",      "Mg++" , "SO4--"}, { -0.0425 } },
    { {"HSO4-",      "Na+"  , "SO4--"}, { -0.0094 } },
    { {"K+",         "Mg++" , "SO4--"}, { -0.0480 } },
    { {"K+",         "Na+"  , "SO4--"}, { -0.0100 } },
    { {"K+",         "OH-"  , "SO4--"}, { -0.0500 } },
    { {"Mg++",       "Na+"  , "SO4--"}, { -0.0150 } },
    { {"Na+",        "OH-"  , "SO4--"}, { -0.0090 } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> lambda_data =
{
    { {"B(OH)3", "Cl-"       }, {  0.091 } },
    { {"B(OH)3", "K+"        }, { -0.140 } },
    { {"B(OH)3", "Na+"       }, { -0.097 } },
    { {"B(OH)3", "SO4--"     }, {  0.018 } },
    { {"B(OH)3", "B3O3(OH)4-"}, { -0.200 } },
    { {"CO2"   , "Ca++"      }, {  0.183 } },
    { {"CO2"   , "Cl-"       }, { -0.005 } },
    { {"CO2"   , "HSO4-"     }, { -0.003 } },
    { {"CO2"   , "K+"        }, {  0.051 } },
    { {"CO2"   , "Mg++"      }, {  0.183 } },
    { {"CO2"   , "Na+"       }, {  0.085 } },
    { {"CO2"   , "SO4--"     }, {  0.097 } },
};

const Pairs<Vec<ChemicalFormula>, Vec<real>> zeta_data =
{
    { {"B(OH)3", "Cl-", "H+"   }, { -0.0102 } },
    { {"B(OH)3", "Na+", "SO4--"}, {  0.0460 } },
    { {"CO2"   , "Na+", "SO4--"}, { -0.0150 } },
};

/// The calculated values of function J0(x) on the points x[i] = 0.0 + i*0.01 in the region (0.0, 1.0)
const Vec<double> J0region1 =
{
    0.00000000e+00, 7.06207558e-05, 2.38544646e-04, 4.80518193e-04, 7.84985974e-04, 1.14430570e-03, 1.55285818e-03, 2.00625822e-03, 2.50094831e-03, 3.03396197e-03,
    3.60277406e-03, 4.20520037e-03, 4.83932706e-03, 5.50345931e-03, 6.19608287e-03, 6.91583456e-03, 7.66147916e-03, 8.43189104e-03, 9.22603924e-03, 1.00429754e-02,
    1.08818235e-02, 1.17417716e-02, 1.26220648e-02, 1.35219989e-02, 1.44409153e-02, 1.53781966e-02, 1.63332630e-02, 1.73055680e-02, 1.82945965e-02, 1.92998614e-02,
    2.03209017e-02, 2.13572800e-02, 2.24085812e-02, 2.34744104e-02, 2.45543916e-02, 2.56481662e-02, 2.67553921e-02, 2.78757422e-02, 2.90089038e-02, 3.01545772e-02,
    3.13124753e-02, 3.24823227e-02, 3.36638550e-02, 3.48568178e-02, 3.60609669e-02, 3.72760670e-02, 3.85018914e-02, 3.97382219e-02, 4.09848476e-02, 4.22415655e-02,
    4.35081791e-02, 4.47844989e-02, 4.60703414e-02, 4.73655292e-02, 4.86698908e-02, 4.99832598e-02, 5.13054751e-02, 5.26363808e-02, 5.39758253e-02, 5.53236619e-02,
    5.66797479e-02, 5.80439451e-02, 5.94161188e-02, 6.07961386e-02, 6.21838773e-02, 6.35792115e-02, 6.49820210e-02, 6.63921888e-02, 6.78096011e-02, 6.92341471e-02,
    7.06657187e-02, 7.21042108e-02, 7.35495206e-02, 7.50015482e-02, 7.64601961e-02, 7.79253690e-02, 7.93969741e-02, 8.08749206e-02, 8.23591201e-02, 8.38494861e-02,
    8.53459342e-02, 8.68483818e-02, 8.83567482e-02, 8.98709545e-02, 9.13909238e-02, 9.29165805e-02, 9.44478508e-02, 9.59846627e-02, 9.75269454e-02, 9.90746298e-02,
    1.00627648e-01, 1.02185934e-01, 1.03749423e-01, 1.05318051e-01, 1.06891756e-01, 1.08470476e-01, 1.10054152e-01, 1.11642726e-01, 1.13236139e-01, 1.14834335e-01, 1.16437258e-01
};

/// The calculated values of function J0(x) on the points x[i] = 1.0 + i*0.1 in the region (1.0, 10.0)
const Vec<double> J0region2 =
{
    1.16437258e-01, 1.32715063e-01, 1.49411766e-01, 1.66486062e-01, 1.83902904e-01, 2.01632240e-01, 2.19648056e-01, 2.37927656e-01, 2.56451090e-01, 2.75200712e-01,
    2.94160824e-01, 3.13317387e-01, 3.32657787e-01, 3.52170640e-01, 3.71845632e-01, 3.91673386e-01, 4.11645345e-01, 4.31753678e-01, 4.51991199e-01, 4.72351292e-01,
    4.92827852e-01, 5.13415235e-01, 5.34108209e-01, 5.54901912e-01, 5.75791820e-01, 5.96773717e-01, 6.17843660e-01, 6.38997965e-01, 6.60233175e-01, 6.81546049e-01,
    7.02933538e-01, 7.24392774e-01, 7.45921053e-01, 7.67515825e-01, 7.89174680e-01, 8.10895338e-01, 8.32675641e-01, 8.54513543e-01, 8.76407104e-01, 8.98354481e-01,
    9.20353920e-01, 9.42403756e-01, 9.64502399e-01, 9.86648338e-01, 1.00884013e+00, 1.03107639e+00, 1.05335582e+00, 1.07567714e+00, 1.09803916e+00, 1.12044072e+00,
    1.14288073e+00, 1.16535811e+00, 1.18787186e+00, 1.21042101e+00, 1.23300461e+00, 1.25562177e+00, 1.27827163e+00, 1.30095335e+00, 1.32366614e+00, 1.34640922e+00,
    1.36918186e+00, 1.39198333e+00, 1.41481295e+00, 1.43767006e+00, 1.46055401e+00, 1.48346418e+00, 1.50639998e+00, 1.52936082e+00, 1.55234615e+00, 1.57535543e+00,
    1.59838813e+00, 1.62144375e+00, 1.64452179e+00, 1.66762179e+00, 1.69074329e+00, 1.71388583e+00, 1.73704899e+00, 1.76023235e+00, 1.78343550e+00, 1.80665805e+00,
    1.82989961e+00, 1.85315981e+00, 1.87643830e+00, 1.89973472e+00, 1.92304873e+00, 1.94638000e+00, 1.96972820e+00, 1.99309303e+00, 2.01647419e+00, 2.03987136e+00, 2.06328427e+00
};

/// The calculated values of function J0(x) on the points x[i] = 10.0 + i*1.0 in the region (10.0, 100.0)
const Vec<double> J0region3 =
{
    2.06328427e+00, 2.29822044e+00, 2.53446209e+00, 2.77181680e+00, 3.01013128e+00, 3.24928140e+00, 3.48916529e+00, 3.72969829e+00, 3.97080934e+00, 4.21243822e+00,
    4.45453343e+00, 4.69705062e+00, 4.93995136e+00, 5.18320209e+00, 5.42677335e+00, 5.67063917e+00, 5.91477653e+00, 6.15916491e+00, 6.40378598e+00, 6.64862328e+00,
    6.89366200e+00, 7.13888874e+00, 7.38429135e+00, 7.62985879e+00, 7.87558100e+00, 8.12144876e+00, 8.36745365e+00, 8.61358792e+00, 8.85984443e+00, 9.10621661e+00,
    9.35269838e+00, 9.59928411e+00, 9.84596857e+00, 1.00927469e+01, 1.03396146e+01, 1.05865675e+01, 1.08336015e+01, 1.10807131e+01, 1.13278988e+01, 1.15751554e+01,
    1.18224798e+01, 1.20698692e+01, 1.23173210e+01, 1.25648325e+01, 1.28124014e+01, 1.30600255e+01, 1.33077027e+01, 1.35554308e+01, 1.38032082e+01, 1.40510328e+01,
    1.42989031e+01, 1.45468174e+01, 1.47947742e+01, 1.50427719e+01, 1.52908093e+01, 1.55388850e+01, 1.57869977e+01, 1.60351462e+01, 1.62833294e+01, 1.65315461e+01,
    1.67797954e+01, 1.70280762e+01, 1.72763875e+01, 1.75247285e+01, 1.77730983e+01, 1.80214960e+01, 1.82699208e+01, 1.85183719e+01, 1.87668486e+01, 1.90153502e+01,
    1.92638760e+01, 1.95124254e+01, 1.97609976e+01, 2.00095922e+01, 2.02582084e+01, 2.05068458e+01, 2.07555038e+01, 2.10041819e+01, 2.12528795e+01, 2.15015962e+01,
    2.17503315e+01, 2.19990850e+01, 2.22478561e+01, 2.24966446e+01, 2.27454499e+01, 2.29942718e+01, 2.32431097e+01, 2.34919634e+01, 2.37408324e+01, 2.39897164e+01, 2.42386152e+01
};

/// The calculated values of function J0(x) on the points x[i] = 100.0 + i*100.0 in the region (100.0, 1000.0)
const Vec<double> J0region4 =
{
    2.42386152e+01, 4.91709893e+01, 7.41388190e+01, 9.91190695e+01, 1.24105388e+02, 1.49095205e+02, 1.74087256e+02, 1.99080833e+02, 2.24075510e+02, 2.49071007e+02
};

/// The calculated values of function J0(x) on the points x[i] = 1000.0 + i*1000.0 in the region (1000.0, 10000.0)
const Vec<double> J0region5 =
{
    2.49071007e+02, 4.99046818e+02, 7.49036364e+02, 9.99030282e+02, 1.24902622e+03, 1.49902328e+03, 1.74902104e+03, 1.99901926e+03, 2.24901780e+03, 2.49901659e+03
};

/// The calculated values of function J1(x) on the points x[i] = 0.0 + i*0.01 in the region (0.0, 1.0)
const Vec<double> J1region1 =
{
    0.00000000e-00, 1.25151744e-04, 4.14749318e-04, 8.24176731e-04, 1.33185071e-03, 1.92376087e-03, 2.58976962e-03, 3.32207863e-03, 4.11444419e-03, 4.96172312e-03,
    5.85958733e-03, 6.80433336e-03, 7.79274934e-03, 8.82201875e-03, 9.88964869e-03, 1.09934151e-02, 1.21313201e-02, 1.33015585e-02, 1.45024900e-02, 1.57326178e-02,
    1.69905698e-02, 1.82750838e-02, 1.95849943e-02, 2.09192222e-02, 2.22767652e-02, 2.36566899e-02, 2.50581251e-02, 2.64802559e-02, 2.79223180e-02, 2.93835939e-02,
    3.08634081e-02, 3.23611240e-02, 3.38761406e-02, 3.54078895e-02, 3.69558327e-02, 3.85194600e-02, 4.00982869e-02, 4.16918532e-02, 4.32997208e-02, 4.49214724e-02,
    4.65567101e-02, 4.82050541e-02, 4.98661416e-02, 5.15396256e-02, 5.32251741e-02, 5.49224690e-02, 5.66312055e-02, 5.83510910e-02, 6.00818448e-02, 6.18231971e-02,
    6.35748887e-02, 6.53366701e-02, 6.71083013e-02, 6.88895510e-02, 7.06801965e-02, 7.24800229e-02, 7.42888229e-02, 7.61063965e-02, 7.79325506e-02, 7.97670984e-02,
    8.16098597e-02, 8.34606598e-02, 8.53193299e-02, 8.71857068e-02, 8.90596322e-02, 9.09409528e-02, 9.28295200e-02, 9.47251900e-02, 9.66278229e-02, 9.85372835e-02,
    1.00453440e-01, 1.02376165e-01, 1.04305334e-01, 1.06240827e-01, 1.08182528e-01, 1.10130321e-01, 1.12084096e-01, 1.14043746e-01, 1.16009167e-01, 1.17980255e-01,
    1.19956912e-01, 1.21939042e-01, 1.23926549e-01, 1.25919343e-01, 1.27917334e-01, 1.29920434e-01, 1.31928560e-01, 1.33941627e-01, 1.35959556e-01, 1.37982267e-01,
    1.40009683e-01, 1.42041730e-01, 1.44078334e-01, 1.46119422e-01, 1.48164927e-01, 1.50214778e-01, 1.52268909e-01, 1.54327256e-01, 1.56389753e-01, 1.58456339e-01, 1.60526953e-01
};

/// The calculated values of function J1(x) on the points x[i] = 1.0 + i*0.1 in the region (1.0, 10.0)
const Vec<double> J1region2 =
{
    1.60526953e-01, 1.81442105e-01, 2.02701762e-01, 2.24262453e-01, 2.46088272e-01, 2.68149211e-01, 2.90419926e-01, 3.12878818e-01, 3.35507329e-01, 3.58289398e-01,
    3.81211036e-01, 4.04259985e-01, 4.27425447e-01, 4.50697864e-01, 4.74068736e-01, 4.97530475e-01, 5.21076279e-01, 5.44700030e-01, 5.68396207e-01, 5.92159812e-01,
    6.15986308e-01, 6.39871564e-01, 6.63811811e-01, 6.87803601e-01, 7.11843772e-01, 7.35929420e-01, 7.60057870e-01, 7.84226655e-01, 8.08433496e-01, 8.32676281e-01,
    8.56953056e-01, 8.81262001e-01, 9.05601427e-01, 9.29969760e-01, 9.54365532e-01, 9.78787370e-01, 1.00323399e+00, 1.02770420e+00, 1.05219687e+00, 1.07671094e+00,
    1.10124543e+00, 1.12579939e+00, 1.15037195e+00, 1.17496229e+00, 1.19956963e+00, 1.22419322e+00, 1.24883238e+00, 1.27348644e+00, 1.29815479e+00, 1.32283683e+00,
    1.34753201e+00, 1.37223979e+00, 1.39695969e+00, 1.42169120e+00, 1.44643390e+00, 1.47118733e+00, 1.49595110e+00, 1.52072481e+00, 1.54550809e+00, 1.57030058e+00,
    1.59510196e+00, 1.61991188e+00, 1.64473005e+00, 1.66955617e+00, 1.69438996e+00, 1.71923115e+00, 1.74407948e+00, 1.76893471e+00, 1.79379660e+00, 1.81866492e+00,
    1.84353946e+00, 1.86842001e+00, 1.89330637e+00, 1.91819835e+00, 1.94309576e+00, 1.96799843e+00, 1.99290619e+00, 2.01781888e+00, 2.04273633e+00, 2.06765841e+00,
    2.09258496e+00, 2.11751584e+00, 2.14245093e+00, 2.16739009e+00, 2.19233319e+00, 2.21728012e+00, 2.24223076e+00, 2.26718500e+00, 2.29214274e+00, 2.31710386e+00, 2.34206827e+00
};

/// The calculated values of function J1(x) on the points x[i] = 10.0 + i*1.0 in the region (10.0, 100.0)
const Vec<double> J1region3 =
{
    2.34206827e+00, 2.59187348e+00, 2.84191919e+00, 3.09214678e+00, 3.34251316e+00, 3.59298607e+00, 3.84354101e+00, 4.09415916e+00, 4.34482586e+00, 4.59552961e+00,
    4.84626133e+00, 5.09701375e+00, 5.34778108e+00, 5.59855865e+00, 5.84934266e+00, 6.10013007e+00, 6.35091839e+00, 6.60170559e+00, 6.85249003e+00, 7.10327038e+00,
    7.35404555e+00, 7.60481466e+00, 7.85557702e+00, 8.10633207e+00, 8.35707936e+00, 8.60781855e+00, 8.85854938e+00, 9.10927165e+00, 9.35998525e+00, 9.61069006e+00,
    9.86138606e+00, 1.01120732e+01, 1.03627516e+01, 1.06134211e+01, 1.08640820e+01, 1.11147342e+01, 1.13653778e+01, 1.16160129e+01, 1.18666397e+01, 1.21172583e+01,
    1.23678687e+01, 1.26184711e+01, 1.28690656e+01, 1.31196523e+01, 1.33702314e+01, 1.36208030e+01, 1.38713673e+01, 1.41219243e+01, 1.43724742e+01, 1.46230171e+01,
    1.48735532e+01, 1.51240825e+01, 1.53746053e+01, 1.56251215e+01, 1.58756314e+01, 1.61261350e+01, 1.63766325e+01, 1.66271241e+01, 1.68776097e+01, 1.71280895e+01,
    1.73785636e+01, 1.76290322e+01, 1.78794953e+01, 1.81299530e+01, 1.83804055e+01, 1.86308528e+01, 1.88812950e+01, 1.91317322e+01, 1.93821645e+01, 1.96325920e+01,
    1.98830149e+01, 2.01334330e+01, 2.03838467e+01, 2.06342558e+01, 2.08846606e+01, 2.11350611e+01, 2.13854573e+01, 2.16358494e+01, 2.18862373e+01, 2.21366213e+01,
    2.23870013e+01, 2.26373774e+01, 2.28877497e+01, 2.31381183e+01, 2.33884832e+01, 2.36388444e+01, 2.38892021e+01, 2.41395563e+01, 2.43899070e+01, 2.46402544e+01, 2.48905984e+01
};

/// The calculated values of function J1(x) on the points x[i] = 100.0 + i*100.0 in the region (100.0, 1000.0)
const Vec<double> J1region4 =
{
    2.48905984e+01, 4.99141212e+01, 7.49270387e+01, 9.99355618e+01, 1.24941747e+02, 1.49946507e+02, 1.74950320e+02, 1.99953466e+02, 2.24956118e+02, 2.49958395e+02
};

/// The calculated values of function J1(x) on the points x[i] = 1000.0 + i*1000.0 in the region (1000.0, 10000.0)
const Vec<double> J1region5 =
{
    2.49958395e+02, 4.99971193e+02, 7.49977057e+02, 9.99980578e+02, 1.24998298e+03, 1.49998474e+03, 1.74998611e+03, 1.99998720e+03, 2.24998810e+03, 2.49998886e+03
};

/// The temperature points (in K) for the interpolation of Aphi from 0 to 350 °C.
const Vec<double> Aphi_temperatures =
{
    273.15, 298.15, 323.15, 348.15, 373.15, 398.15, 423.15, 448.15, 473.15, 498.15, 523.15, 548.15, 573.15, 598.15, 623.15
};

/// The pressure points (in Pa) for the interpolation of Aphi from 1 to 5000 bar.
const Vec<double> Aphi_pressures =
{
    1e5, 100e5, 200e5, 400e5, 600e5, 800e5, 1000e5, 1500e5, 2000e5, 3000e5, 4000e5, 5000e5,
};

/// The data used for the interpolation of Aphi.
/// The following equation was used to compute these data:
/// \f[
///      A^{\phi}=1.400684\cdot10^{6}\frac{1}{\rho_{w}}\left(\frac{\rho_{w}}{\epsilon_{w}T}\right)^{\frac{3}{2}},
/// \f]
/// where @f$T@f$ is temperature in K, @f$\epsilon_{w}@f$ is the dielectric constant of water, and @f$\rho_{w}@f$ its density in g/m<sup>3<\sup>
const Vec<double> Aphi_data =
{
    0.37674222, 0.39147278, 0.41031975, 0.43330198, 0.46058910, 0.49248630, 0.52954906, 0.57258689, 0.62280952, 0.68208055, 0.75343848, 0.84226817, 0.95932062, 1.13030749, 1.43881804,
    0.37511283, 0.38962104, 0.40816127, 0.43071699, 0.45741532, 0.48854946, 0.52461006, 0.56635820, 0.61497721, 0.67238823, 0.74195091, 0.83021165, 0.95211387, 1.13030749, 1.43881804,
    0.37349119, 0.38778969, 0.40603633, 0.42818435, 0.45432255, 0.48468653, 0.51967722, 0.55991097, 0.60632611, 0.66040079, 0.72461917, 0.80356970, 0.90695201, 1.06002794, 1.35961218,
    0.37031872, 0.38423960, 0.40194465, 0.42334123, 0.44845489, 0.47742645, 0.51051253, 0.54810417, 0.59077384, 0.63937254, 0.69522720, 0.76056419, 0.83951064, 0.94084255, 1.08752618,
    0.36723706, 0.38083131, 0.39805034, 0.41877218, 0.44297410, 0.47072404, 0.50217061, 0.53754295, 0.57716443, 0.62148709, 0.67116338, 0.72719329, 0.79124059, 0.86637464, 0.95905835,
    0.36424217, 0.37755559, 0.39433845, 0.41445319, 0.43784086, 0.46451346, 0.49453834, 0.52802729, 0.56513171, 0.60604478, 0.65101410, 0.70037330, 0.75461241, 0.81453466, 0.88162611,
    0.36133042, 0.37440405, 0.39079550, 0.41036291, 0.43302090, 0.45873872, 0.48752239, 0.51939798, 0.55439682, 0.59254133, 0.63383034, 0.67822265, 0.72561799, 0.77583430, 0.82858092,
    0.35439502, 0.36702034, 0.38260006, 0.40101786, 0.42215421, 0.44591174, 0.47220102, 0.50092041, 0.53193360, 0.56504392, 0.59996245, 0.63626628, 0.67334062, 0.71029469, 0.74583251,
    0.34791990, 0.36026613, 0.37522638, 0.39274277, 0.41269203, 0.43494673, 0.45937116, 0.48580434, 0.51403915, 0.54379790, 0.57470298, 0.60623985, 0.63770811, 0.66815393, 0.69627074,
    0.33621768, 0.34833639, 0.36246295, 0.37869462, 0.39694660, 0.41708490, 0.43894731, 0.46233761, 0.48701196, 0.51266165, 0.53889233, 0.56519854, 0.59093075, 0.61524994, 0.63706014,
    0.32597289, 0.33811759, 0.35176187, 0.36716097, 0.38429206, 0.40304427, 0.42326335, 0.44475696, 0.46728781, 0.49056178, 0.51421251, 0.53778186, 0.56069387, 0.58221783, 0.60141144,
    0.31695465, 0.32925197, 0.34262585, 0.35747187, 0.37383264, 0.39162945, 0.41072659, 0.43094633, 0.45206745, 0.47381721, 0.49585921, 0.51777702, 0.53905156, 0.55902802, 0.57686399
};

/// Return the coefficients for the Pitzer parameter function of a groups of
/// chemical formulas. These coefficients are searched in given `data`. If
/// `data` does not contain a row with chemical formulas equivalent to the
/// given ones in `group`, an empty vector of coefficients are returned.
/// Otherwise, the corresponding coefficients are returned.
auto getParamFunctionCoeffs(Vec<ChemicalFormula> const& group, Pairs<Vec<ChemicalFormula>, Vec<real>> const& data) -> Index
{
    // Auxiliary function that returns true if two vectors of **sorted chemical formulas** are component-wise equivalent.
    auto equivalent_groups = [](auto const& group1, auto const& group2)
    {
        assert(group1.size() != group2.size());
        for(auto i = 0; i < group1.size(); ++i)
            if(!ChemicalFormula::equivalent(group1[i], group2[i]))
                return false;
        return true;
    };

    auto const sortedgroup = sorted(group);

    for(auto const& [i, line] : enumerate(data))
    {
        auto const& [othergroup, coeffs] = line; // othergroup is expected to have been sorted at construction of data!

        if(equivalent_groups(sortedgroup, othergroup))
            return i;
    }

    return data.size();
}

/// Return a function that computes a Pitzer interaction parameter at a given temperature (in K).
/// @param c The vector with coefficients for the Pitzer interaction parameter function.
auto createInteractionParamFunctionWithCoeffs(Vec<real> const& c) -> Fn<real(real)>
{
    auto const Tr = 298.15;

    if(c.size() == 0) return [=](real T) { return 0.0; };
    if(c.size() == 1) return [=](real T) { return c[0]; };
    if(c.size() == 2) return [=](real T) { return c[0] + c[1]*(1.0/T - 1.0/Tr); };
    if(c.size() == 3) return [=](real T) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr); };
    if(c.size() == 4) return [=](real T) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr); };
    if(c.size() == 5) return [=](real T) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr); };
    if(c.size() == 6) return [=](real T) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr) + c[5]*(1.0/(T*T) - 1.0/(Tr*Tr)); };

    errorif(true, "Cannot create the temperature-dependent parameter function with given coefficients: ", str(c));
}

/// Return a function that computes a Pitzer interaction parameter at a given temperature (in K) for a given species group.
/// @param group The list of chemical formulas specifying the species group.
/// @param data The table containing groups of species formulas and their respective coefficients for their interaction parameter function.
auto createInteractionParamFunctionForSpeciesGroup(Vec<ChemicalFormula> const& group, Pairs<Vec<ChemicalFormula>, Vec<real>> const& data) -> Fn<real(real)>
{
    // Auxiliary function that returns true if two vectors of **sorted chemical formulas** are component-wise equivalent.
    auto equivalent_groups = [](auto const& group1, auto const& group2)
    {
        assert(group1.size() != group2.size());
        for(auto i = 0; i < group1.size(); ++i)
            if(!ChemicalFormula::equivalent(group1[i], group2[i]))
                return false;
        return true;
    };

    auto const sortedgroup = sorted(group);

    for(auto const& line : data)
    {
        auto const& [othergroup, coeffs] = line; // othergroup is expected to have been sorted at construction of data!
        if(equivalent_groups(sortedgroup, othergroup))
            return createInteractionParamFunctionWithCoeffs(coeffs);
    }

    return createInteractionParamFunctionWithCoeffs({});
}

auto createTableInteractionParamFunctionsBinary(Vec<ChemicalFormula> const& formulas1, Vec<ChemicalFormula> const& formulas2, Pairs<Vec<ChemicalFormula>, Vec<real>> const& data) -> Table2D<Fn<real(real)>>
{
    Table2D<Fn<real(real)>> res = initTable2D<Fn<real(real)>>(formulas1.size(), formulas2.size());

    for(auto i = 0; i < formulas1.size(); ++i)
        for(auto j = 0; j < formulas2.size(); ++j)
            res[i][j] = createInteractionParamFunctionForSpeciesGroup({formulas1[i], formulas2[j]}, data);

    return res;
}

auto createTableInteractionParamFunctionsTernary(Vec<ChemicalFormula> const& formulas1, Vec<ChemicalFormula> const& formulas2, Vec<ChemicalFormula> const& formulas3, Pairs<Vec<ChemicalFormula>, Vec<real>> const& data) -> Table3D<Fn<real(real)>>
{
    Table3D<Fn<real(real)>> res = initTable3D<Fn<real(real)>>(formulas1.size(), formulas2.size(), formulas3.size());

    for(auto i = 0; i < formulas1.size(); ++i)
        for(auto j = 0; j < formulas2.size(); ++j)
            for(auto k = 0; k < formulas3.size(); ++j)
                res[i][j][k] = createInteractionParamFunctionForSpeciesGroup({formulas1[i], formulas2[j], formulas3[k]}, data);

    return res;
}

auto createBeta0Table(Vec<ChemicalFormula> const& cations, Vec<ChemicalFormula> const& anions) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(cations, anions, beta0_data);
}

auto createBeta1Table(Vec<ChemicalFormula> const& cations, Vec<ChemicalFormula> const& anions) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(cations, anions, beta1_data);
}

auto createBeta2Table(Vec<ChemicalFormula> const& cations, Vec<ChemicalFormula> const& anions) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(cations, anions, beta2_data);
}

auto createCphiTable(Vec<ChemicalFormula> const& cations, Vec<ChemicalFormula> const& anions) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(cations, anions, Cphi_data);
}

auto createThetaTable(Vec<ChemicalFormula> const& ions1, Vec<ChemicalFormula> const& ions2) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(ions1, ions2, theta_data);
}

auto createPsiTable(Vec<ChemicalFormula> const& ions1, Vec<ChemicalFormula> const& ions2, Vec<ChemicalFormula> const& ions3) -> Table3D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsTernary(ions1, ions2, ions3, psi_data);
}

auto createLambdaTable(Vec<ChemicalFormula> const& neutrals, Vec<ChemicalFormula> const& ions) -> Table2D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsBinary(neutrals, ions, lambda_data);
}

auto createZetaTable(Vec<ChemicalFormula> const& neutrals, Vec<ChemicalFormula> const& cations, Vec<ChemicalFormula> const& anions) -> Table3D<Fn<real(real)>>
{
    return createTableInteractionParamFunctionsTernary(neutrals, cations, anions, zeta_data);
}

auto interpolate(real x, double x0, double x1, Vec<double> const& ypoints) -> real
{
    const auto N  = ypoints.size();
    const auto dx = (x1 - x0)/(N - 1);

    const auto i = Index((x - x0)/dx);

    return (i < 0) ? ypoints.front() : (i > N) ? ypoints.back() : ypoints[i];
}

auto J0(real x) -> real
{
    if(0.0 <= x && x <= 1.0) return interpolate(x, 0.0, 1.0, J0region1);
    if(1.0 <= x && x <= 10.0) return interpolate(x, 1.0, 10.0, J0region2);
    if(10.0 <= x && x <= 100.0) return interpolate(x, 10.0, 100.0, J0region3);
    if(100.0 <= x && x <= 1000.0) return interpolate(x, 100.0, 1000.0, J0region4);
    if(1000.0 <= x && x <= 10000.0) return interpolate(x, 1000.0, 10000.0, J0region5);

    Exception exception;
    exception.error << "Cannot interpolate the Pitzer function J0(x) with the provided x = " << x << ".";
    exception.reason << "The value of x must be within the interval [0, 10000].";
    RaiseError(exception);

    return 0.0;
}

auto J1(real x) -> real
{
    if(0.0 <= x && x <= 1.0) return interpolate(x, 0.0, 1.0, J1region1);
    if(1.0 <= x && x <= 10.0) return interpolate(x, 1.0, 10.0, J1region2);
    if(10.0 <= x && x <= 100.0) return interpolate(x, 10.0, 100.0, J1region3);
    if(100.0 <= x && x <= 1000.0) return interpolate(x, 100.0, 1000.0, J1region4);
    if(1000.0 <= x && x <= 10000.0) return interpolate(x, 1000.0, 10000.0, J1region5);

    Exception exception;
    exception.error << "Cannot interpolate the Pitzer function J1(x) with the provided x = " << x << ".";
    exception.reason << "The value of x must be within the interval [0, 10000].";
    RaiseError(exception);

    return 0.0;
}

struct PitzerParamsState
{
    Table2D<real> beta0;
    Table2D<real> beta1;
    Table2D<real> beta2;
    Table2D<real> Cphi;
    Table2D<real> theta_cc;
    Table2D<real> theta_aa;
    Table3D<real> psi_cca;
    Table3D<real> psi_aac;
    Table2D<real> lambda_nc;
    Table2D<real> lambda_na;
    Table3D<real> zeta;
    real Aphi;

    PitzerParamsState() = default;

    explicit PitzerParamsState(Index Nc, Index Na, Index Nn)
    {
        beta0     = initTable2D<real>(Nc, Na);
        beta1     = initTable2D<real>(Nc, Na);
        beta2     = initTable2D<real>(Nc, Na);
        Cphi      = initTable2D<real>(Nc, Na);
        theta_cc  = initTable2D<real>(Nc, Nc);
        theta_aa  = initTable2D<real>(Na, Na);
        psi_cca   = initTable3D<real>(Nc, Nc, Na);
        psi_aac   = initTable3D<real>(Na, Na, Nc);
        lambda_nc = initTable2D<real>(Nn, Nc);
        lambda_na = initTable2D<real>(Nn, Na);
        zeta      = initTable3D<real>(Nn, Nc, Na);
    }
};

struct PitzerParams
{
    PitzerParams() = default;

    explicit PitzerParams(AqueousMixture const& solution);

    const Vec<ChemicalFormula> neutrals;     ///< The chemical formulas of the neutral species in the aqueous solution.
    const Vec<ChemicalFormula> charged;      ///< The chemical formulas of the charged species in the aqueous solution.
    const Vec<ChemicalFormula> cations;      ///< The chemical formulas of the cations in the aqueous solution.
    const Vec<ChemicalFormula> anions;       ///< The chemical formulas of the anions in the aqueous solution.
    const Index Nn;                          ///< The number of neutral species in the aqueous solution.
    const Index Nz;                          ///< The number of charged species in the aqueous solution.
    const Index Nc;                          ///< The number of cations in the aqueous solution.
    const Index Na;                          ///< The number of anions in the aqueous solution.
    const Indices idx_neutrals;              ///< The indices of the neutral species among all species in the aqueous solution.
    const Indices idx_charged;               ///< The indices of the charged species among all species in the aqueous solution.
    const Indices idx_cations;               ///< The indices of the cations among all species in the aqueous solution.
    const Indices idx_anions;                ///< The indices of the anions among all species in the aqueous solution.
    const ArrayXr z_charged;                 ///< The electric charge of the charged species.
    const ArrayXr z_cations;                 ///< The electric charge of the cations.
    const ArrayXr z_anions;                  ///< The electric charge of the anions.
    const Table2D<Fn<real(real)>> beta0;     ///< The functions \eq{\beta^{(0)}_{MX}(T, P)} in the Pitzer model.
    const Table2D<Fn<real(real)>> beta1;     ///< The functions \eq{\beta^{(1)}_{MX}(T, P)} in the Pitzer model.
    const Table2D<Fn<real(real)>> beta2;     ///< The functions \eq{\beta^{(2)}_{MX}(T, P)} in the Pitzer model.
    const Table2D<Fn<real(real)>> Cphi;      ///< The functions \eq{C^{\phi}_{MX}(T, P)} in the Pitzer model.
    const Table2D<Fn<real(real)>> theta_cc;  ///< The functions \eq{\theta_{ij}(T, P)} in the Pitzer model for cation-cation interactions.
    const Table2D<Fn<real(real)>> theta_aa;  ///< The functions \eq{\theta_{ij}(T, P)} in the Pitzer model for anion-anion interactions.
    const Table3D<Fn<real(real)>> psi_cca;   ///< The functions \eq{\psi_{ijk}(T, P)} in the Pitzer model for cation-cation-anion interactions.
    const Table3D<Fn<real(real)>> psi_aac;   ///< The functions \eq{A^\phi(T, P)} in the Pitzer model for anion-anion-cation interactions.
    const Table2D<Fn<real(real)>> lambda_nc; ///< The functions \eq{\lambda_{ij}(T, P)} in the Pitzer model for neutral-cation interactions.
    const Table2D<Fn<real(real)>> lambda_na; ///< The functions \eq{\lambda_{ij}(T, P)} in the Pitzer model for neutral-anion interactions.
    const Table3D<Fn<real(real)>> zeta;      ///< The functions \eq{\zeta^\phi(T, P)} in the Pitzer model for neutral-cation-anion interactions.
    const BilinearInterpolator Aphi;         ///< The function \eq{A^\phi(T, P)} in the Pitzer model.

    PitzerParamsState state;
};

PitzerParams::PitzerParams(AqueousMixture const& solution)
: neutrals( vectorize(solution.neutral(), RKT_LAMBDA(x, x.formula())) ),
  charged( vectorize(solution.charged(), RKT_LAMBDA(x, x.formula())) ),
  cations( vectorize(solution.cations(), RKT_LAMBDA(x, x.formula())) ),
  anions( vectorize(solution.anions(), RKT_LAMBDA(x, x.formula())) ),
  Nn( neutrals.size() ),
  Nz( charged.size() ),
  Nc( cations.size() ),
  Na( anions.size() ),
  idx_neutrals( solution.indicesNeutral() ),
  idx_charged( solution.indicesCharged() ),
  idx_cations( solution.indicesCations() ),
  idx_anions( solution.indicesAnions() ),
  z_charged( solution.charges()(idx_charged) ),
  z_cations( solution.charges()(idx_cations) ),
  z_anions( solution.charges()(idx_anions) ),
  beta0( createBeta0Table(cations, anions) ),
  beta1( createBeta1Table(cations, anions) ),
  beta2( createBeta2Table(cations, anions) ),
  Cphi ( createCphiTable(cations, anions) ),
  theta_cc( createThetaTable(cations, cations) ),
  theta_aa( createThetaTable(anions, anions) ),
  psi_cca( createPsiTable(cations, cations, anions) ),
  psi_aac( createPsiTable(anions, anions, cations) ),
  lambda_nc( createLambdaTable(neutrals, cations) ),
  lambda_na( createLambdaTable(neutrals, anions) ),
  zeta( createZetaTable(neutrals, cations, anions) ),
  Aphi( BilinearInterpolator(Aphi_temperatures, Aphi_pressures, Aphi_data) ),
  state(Nc, Na, Nn)
{
}

auto thetaE(AqueousMixtureState const& state, PitzerParams const& pitzer, real zi, real zj) -> real
{
    if(zi == zj) return 0.0;

    const auto I     = state.Ie;
    const auto T     = state.T;
    const auto P     = state.P;
    const auto sqrtI = sqrt(I);
    const auto Aphi  = pitzer.state.Aphi;
    const auto xij   = 6.0*zi*zj*Aphi*sqrtI;
    const auto xii   = 6.0*zi*zi*Aphi*sqrtI;
    const auto xjj   = 6.0*zj*zj*Aphi*sqrtI;
    const auto J0ij  = J0(xij);
    const auto J0ii  = J0(xii);
    const auto J0jj  = J0(xjj);

    return zi*zj/(4*I) * (J0ij - 0.5*J0ii - 0.5*J0jj);
}

auto thetaE_prime(AqueousMixtureState const& state, PitzerParams const& pitzer, real zi, real zj) -> real
{
    if(zi == zj) return 0.0;

    const auto I     = state.Ie;
    const auto T     = state.T;
    const auto P     = state.P;
    const auto sqrtI = sqrt(I);
    const auto Aphi  = pitzer.state.Aphi;
    const auto xij   = 6.0*zi*zj*Aphi*sqrtI;
    const auto xii   = 6.0*zi*zi*Aphi*sqrtI;
    const auto xjj   = 6.0*zj*zj*Aphi*sqrtI;
    const auto J1ij  = J1(xij);
    const auto J1ii  = J1(xii);
    const auto J1jj  = J1(xjj);

    return zi*zj/(8*I*I) * (J1ij - 0.5*J1ii - 0.5*J1jj) - thetaE(state, pitzer, zi, zj)/I;
}

auto thetaE_cc(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto zi = pitzer.z_cations[i];
    const auto zj = pitzer.z_cations[j];
    return thetaE(state, pitzer, zi, zj);
}

auto thetaE_aa(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto zi = pitzer.z_anions[i];
    const auto zj = pitzer.z_anions[j];
    return thetaE(state, pitzer, zi, zj);
}

auto thetaE_prime_cc(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto zi = pitzer.z_cations[i];
    const auto zj = pitzer.z_cations[j];
    return thetaE_prime(state, pitzer, zi, zj);
}

auto thetaE_prime_aa(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto zi = pitzer.z_anions[i];
    const auto zj = pitzer.z_anions[j];
    return thetaE_prime(state, pitzer, zi, zj);
}

auto Phi_cc(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto thetaij = pitzer.state.theta_cc[i][j];
    return thetaij + thetaE_cc(state, pitzer, i, j);
}

auto Phi_aa(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto thetaij = pitzer.state.theta_aa[i][j];
    return thetaij + thetaE_aa(state, pitzer, i, j);
}

auto Phi_phi_cc(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto I = state.Ie;
    const auto thetaij = pitzer.state.theta_cc[i][j];
    return thetaij + thetaE_cc(state, pitzer, i, j) + I * thetaE_prime_cc(state, pitzer, i, j);
}

auto Phi_phi_aa(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    const auto I = state.Ie;
    const auto thetaij = pitzer.state.theta_aa[i][j];
    return thetaij + thetaE_aa(state, pitzer, i, j) + I * thetaE_prime_aa(state, pitzer, i, j);
}

auto Phi_prime_cc(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    return thetaE_prime_cc(state, pitzer, i, j);
}

auto Phi_prime_aa(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned i, unsigned j) -> real
{
    return thetaE_prime_aa(state, pitzer, i, j);
}

auto g(real x) -> real
{
    return 2.0*(1 - (1 + x)*exp(-x))/(x*x);
}

auto g_prime(real x) -> real
{
    return -2.0*(1 - (1 + x + 0.5*x*x)*exp(-x))/(x*x);
}

const auto alpha  =  2.0;
const auto alpha1 =  1.4;
const auto alpha2 = 12.0;

auto B_phi(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned c, unsigned a) -> real
{
    const auto I     = state.Ie;
    const auto T     = state.T;
    const auto sqrtI = sqrt(I);
    const auto zc    = pitzer.z_cations[c];
    const auto za    = pitzer.z_anions[a];
    const auto beta0 = pitzer.beta0[c][a](T);
    const auto beta1 = pitzer.beta1[c][a](T);
    const auto beta2 = pitzer.beta2[c][a](T);

    if(abs(zc) == 2 && abs(za) == 2)
        return beta0 + beta1 * exp(-alpha1*sqrtI) + beta2 * exp(-alpha2*sqrtI);
    else
        return beta0 + beta1 * exp(-alpha*sqrtI);
}

auto B(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned c, unsigned a) -> real
{
    const auto I     = state.Ie;
    const auto T     = state.T;
    const auto sqrtI = sqrt(I);
    const auto zc    = pitzer.z_cations[c];
    const auto za    = pitzer.z_anions[a];
    const auto beta0 = pitzer.beta0[c][a](T);
    const auto beta1 = pitzer.beta1[c][a](T);
    const auto beta2 = pitzer.beta2[c][a](T);

    if(abs(zc) == 2 && abs(za) == 2)
        return beta0 + beta1 * g(alpha1*sqrtI) + beta2 * g(alpha2*sqrtI);
    else
        return beta0 + beta1 * g(alpha*sqrtI);
}

auto B_prime(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned c, unsigned a) -> real
{
    const auto I     = state.Ie;
    const auto T     = state.T;
    const auto sqrtI = sqrt(I);
    const auto zc    = pitzer.z_cations[c];
    const auto za    = pitzer.z_anions[a];
    const auto beta1 = pitzer.beta1[c][a](T);
    const auto beta2 = pitzer.beta2[c][a](T);

    if(abs(zc) == 2 && abs(za) == 2)
        return beta1 * g_prime(alpha1*sqrtI)/I + beta2 * g_prime(alpha2*sqrtI)/I;
    else
        return beta1 * g_prime(alpha*sqrtI)/I;
}

auto C(AqueousMixtureState const& state, PitzerParams const& pitzer, unsigned c, unsigned a) -> real
{
    const auto T     = state.T;
    const auto Cphi = pitzer.Cphi[c][a](T);
    const auto zc   = pitzer.z_cations[c];
    const auto za   = pitzer.z_anions[a];
    return 0.5 * Cphi/sqrt(abs(zc*za));
}

auto computeF(AqueousMixtureState const& state, PitzerParams const& pitzer) -> real
{
    // The indices of the cations and anions
    auto const& idx_cations = pitzer.idx_cations;
    auto const& idx_anions  = pitzer.idx_anions;

    // The number of the cations and anions
    const auto num_cations = idx_cations.size();
    const auto num_anions  = idx_anions.size();

    // THe temperature and pressure in units of K and Pa respectively
    auto const& T = state.T;
    auto const& P = state.P;

    // The molalities of all aqueous species
    auto const& m = state.m;

    // The ionic strength of the aqueous solution and its square root
    const auto I = state.Ie;
    const auto sqrtI = sqrt(I);

    // The Debye-Huckel coefficient Aphi
    const auto Aphi = pitzer.state.Aphi;

    // The b parameter of the Harvie-Moller-Weare Pitzer's model
    const auto b = 1.2;

    // Calculate the term F of the Harvie-Moller-Weare Pitzer's model
    real F = -Aphi * (sqrtI/(1 + b*sqrtI) + 2.0/b * log(1 + b*sqrtI));

    // Iterate over all pairs of cations and anions
    for(auto c = 0; c < num_cations; ++c) for(auto a = 0; a < num_anions; ++a)
        F += m[idx_cations[c]] * m[idx_anions[a]] * B_prime(state, pitzer, c, a);

    // Iterate over all pairs of distinct cations
    for(auto i = 0; i < num_cations - 1; ++i) for(auto j = i + 1; j < num_cations; ++j)
        F += m[idx_cations[i]] * m[idx_cations[j]] * Phi_prime_cc(state, pitzer, i, j);

    // Iterate over all pairs of distinct anions
    for(auto i = 0; i < num_anions - 1; ++i) for(auto j = i + 1; j < num_anions; ++j)
        F += m[idx_anions[i]] * m[idx_anions[j]] * Phi_prime_aa(state, pitzer, i, j);

    return F;
}

auto computeZ(AqueousMixtureState const& state, PitzerParams const& pitzer) -> real
{
    const auto mi = state.m(pitzer.idx_charged);
    const auto zi = pitzer.z_charged.abs();
    return (mi * zi).sum();
}

/// Return the Pitzer activity coefficient of a cation (in natural log scale).
/// @param state The state of the aqueous solution
/// @param pitzer The Pitzer parameters
/// @param M The local index of the cation among all cations in the solution
auto lnActivityCoefficientCation(AqueousMixtureState const& state, PitzerParams const& pitzer, Index M) -> real
{
    // The indices of the neutral species, cations and anions
    auto const& idx_neutrals = pitzer.idx_neutrals;
    auto const& idx_cations  = pitzer.idx_cations;
    auto const& idx_anions   = pitzer.idx_anions;

    // The number of neutrals, cations and anions
    const auto num_neutrals = idx_neutrals.size();
    const auto num_cations  = idx_cations.size();
    const auto num_anions   = idx_anions.size();

    // The vector of molalities of all aqueous species
    auto const& m = state.m;

    // The number of species in the solution
    const auto nspecies = m.size();

    // The electrical charge of the M-th cation
    const auto zM = pitzer.z_cations[M];

    // The terms F and Z of the Harvie-Moller-Weare Pitzer's model
    const auto F = computeF(state, pitzer);
    const auto Z = computeZ(state, pitzer);

    // The log of the activity coefficient of the M-th cation
    real ln_gammaM = {};

    // Iterate over all anions
    for(auto a = 0; a < num_anions; ++a)
    {
        const auto ma = m[idx_anions[a]];
        const auto BMa = B(state, pitzer, M, a);
        const auto CMa = C(state, pitzer, M, a);

        ln_gammaM += ma * (2*BMa + Z*CMa);
    }

    // Iterate over all cations
    for(auto c = 0; c < num_cations; ++c)
    {
        const auto mc = m[idx_cations[c]];
        const auto PhiMc = Phi_cc(state, pitzer, M, c);

        real aux = 2*PhiMc;

        // Iterate over all anions
        for(auto a = 0; a < num_anions; ++a)
        {
            const auto ma = m[idx_anions[a]];
            const auto psiMca = pitzer.state.psi_cca[M][c][a];

            aux += ma * psiMca;
        }

        ln_gammaM += mc * aux;
    }

    // Iterate over all pairs of distinct anions
    for(auto i = 0; i < num_anions - 1; ++i) for(auto j = i + 1; j < num_anions; ++j)
    {
        const auto mi = m[idx_anions[i]];
        const auto mj = m[idx_anions[j]];
        const auto psi_ijM = pitzer.state.psi_aac[i][j][M];

        ln_gammaM += mi * mj * psi_ijM;
    }

    // Iterate over all pairs of cations and anions
    for(auto c = 0; c < num_cations; ++c) for(auto a = 0; a < num_anions; ++a)
    {
        const auto mc = m[idx_cations[c]];
        const auto ma = m[idx_anions[a]];
        const auto Cca = C(state, pitzer, c, a);

        ln_gammaM += abs(zM) * mc * ma * Cca;
    }

    // Iterate over all pairs of neutral species
    for(auto n = 0; n < num_neutrals; ++n)
        ln_gammaM += 2.0 * m[idx_neutrals[n]] * pitzer.state.lambda_nc[n][M];

    // Finalize the calculation
    ln_gammaM += zM*zM*F;

    return ln_gammaM;
}

/// Return the Pitzer activity coefficient of an anion (in natural log scale).
/// @param state The state of the aqueous solution
/// @param pitzer The Pitzer parameters
/// @param X The local index of the anion among all anions in the solution
auto lnActivityCoefficientAnion(AqueousMixtureState const& state, PitzerParams const& pitzer, Index X) -> real
{
    // The indices of the neutral species, cations and anions
    auto const& idx_neutrals = pitzer.idx_neutrals;
    auto const& idx_cations  = pitzer.idx_cations;
    auto const& idx_anions   = pitzer.idx_anions;

    // The number of neutrals, cations and anions
    const auto num_neutrals = idx_neutrals.size();
    const auto num_cations  = idx_cations.size();
    const auto num_anions   = idx_anions.size();

    // The molalities of all aqueous species
    auto const& m = state.m;

    // The number of species in the solution
    const auto nspecies = m.size();

    // The electrical charge of the X-th anion
    const auto zX = pitzer.z_anions[X];

    // The terms F and Z of the Harvie-Moller-Weare Pitzer's model
    const auto F = computeF(state, pitzer);
    const auto Z = computeZ(state, pitzer);

    // The log of the activity coefficient of the X-th anion
    real ln_gammaX = {};

    // Iterate over all cations
    for(auto c = 0; c < num_cations; ++c)
    {
        const auto mc = m[idx_cations[c]];
        const auto BcX = B(state, pitzer, c, X);
        const auto CcX = C(state, pitzer, c, X);

        ln_gammaX += mc * (2*BcX + Z*CcX);
    }

    // Iterate over all anions
    for(auto a = 0; a < num_anions; ++a)
    {
        const auto ma = m[idx_anions[a]];
        const auto PhiXa = Phi_aa(state, pitzer, X, a);

        real aux = 2*PhiXa;

        // Iterate over all cations
        for(auto c = 0; c < num_cations; ++c)
        {
            const auto mc = m[idx_cations[c]];
            const auto psiXac = pitzer.state.psi_aac[X][a][c];

            aux += mc * psiXac;
        }

        ln_gammaX += ma * aux;
    }

    // Iterate over all pairs of distinct cations
    for(auto i = 0; i < num_cations - 1; ++i) for(auto j = i + 1; j < num_cations; ++j)
    {
        const auto mi = m[idx_cations[i]];
        const auto mj = m[idx_cations[j]];
        const auto psi_ijX = pitzer.state.psi_cca[i][j][X];

        ln_gammaX += mi * mj * psi_ijX;
    }

    // Iterate over all pairs of cations and anions
    for(auto c = 0; c < num_cations; ++c) for(auto a = 0; a < num_anions; ++a)
    {
        const auto mc = m[idx_cations[c]];
        const auto ma = m[idx_anions[a]];
        const auto Cca = C(state, pitzer, c, a);

        ln_gammaX += abs(zX) * mc * ma * Cca;
    }

    // Iterate over all pairs of neutral species
    for(auto n = 0; n < num_neutrals; ++n)
        ln_gammaX += 2.0 * m[idx_neutrals[n]] * pitzer.state.lambda_na[n][X];

    // Finalize the calculation
    ln_gammaX += zX*zX*F;

    return ln_gammaX;
}

/// Return the Pitzer activity of water (in natural log scale).
/// @param state The state of the aqueous solution
/// @param pitzer The Pitzer parameters
/// @param iH2O The index of the water species
auto lnActivityWater(AqueousMixtureState const& state, PitzerParams const& pitzer, Index iH2O) -> real
{
    // The indices of the neutral species, cations and anions
    auto const& idx_neutrals = pitzer.idx_neutrals;
    auto const& idx_cations  = pitzer.idx_cations;
    auto const& idx_anions   = pitzer.idx_anions;

    // The number of neutrals, cations and anions
    const auto num_neutrals = idx_neutrals.size();
    const auto num_cations  = idx_cations.size();
    const auto num_anions   = idx_anions.size();

    // The vector of molalities of all aqueous species
    auto const& m = state.m;

    // Calculate the sum of molalities of the solutes
    const auto sum_mi = sum(m) - m[iH2O];

    // Check if pure water conditions (or ions in tiny amounts compared to the amount of water) to avoid division by zero below!
    if(sum_mi == 0.0)
        return 0.0;

    // The number of species in the solution
    const auto nspecies = m.size();

    // The ionic strength of the aqueous solution
    const auto I = state.Ie;

    // The square root of the ionic strength of the aqueous solution
    const auto sqrtI = sqrt(I);

    // The molar mass of water
    const auto Mw = waterMolarMass;

    // The Debye-Huckel coefficient Aphi
    const auto Aphi = pitzer.state.Aphi;

    // The b parameter of the Harvie-Moller-Weare Pitzer's model
    const auto b = 1.2;

    // The term Z of the Harvie-Moller-Weare Pitzer's model
    const auto Z = computeZ(state, pitzer);

    // The osmotic coefficient of the aqueous solution
    real phi = -Aphi*I*sqrtI/(1 + b*sqrtI);

    // Iterate over all pairs of cations and anions
    for(auto c = 0; c < num_cations; ++c) for(auto a = 0; a < num_anions; ++a)
    {
        const auto mc = m[idx_cations[c]];
        const auto ma = m[idx_anions[a]];
        const auto Bca = B_phi(state, pitzer, c, a);
        const auto Cca = C(state, pitzer, c, a);

        phi += mc * ma * (Bca + Z*Cca);
    }

    // Iterate over all pairs of distinct cations
    for(auto i = 0; i < num_cations - 1; ++i) for(auto j = i + 1; j < num_cations; ++j)
    {
        const auto mi = m[idx_cations[i]];
        const auto mj = m[idx_cations[j]];
        const auto Phi_ij = Phi_phi_cc(state, pitzer, i, j);

        real aux = Phi_ij;

        // Iterate over all anions
        for(auto a = 0; a < num_anions; ++a)
        {
            const auto ma = m[idx_anions[a]];
            const auto psi_ija = pitzer.state.psi_cca[i][j][a];

            aux += ma * psi_ija;
        }

        phi += mi * mj * aux;
    }

    // Iterate over all pairs of distinct anions
    for(auto i = 0; i < num_anions - 1; ++i) for(auto j = i + 1; j < num_anions; ++j)
    {
        const auto mi = m[idx_anions[i]];
        const auto mj = m[idx_anions[j]];
        const auto Phi_ij = Phi_phi_aa(state, pitzer, i, j);

        real aux = Phi_ij;

        // Iterate over all cations
        for(auto c = 0; c < num_cations; ++c)
        {
            const auto mc = m[idx_cations[c]];
            const auto psi_ija = pitzer.state.psi_aac[i][j][c];

            aux += mc * psi_ija;
        }

        phi += mi * mj * aux;
    }

    // Iterate over all pairs of neutral species and cations
    for(auto n = 0; n < num_neutrals; ++n)
        for(auto c = 0; c < num_cations; ++c)
            phi += m[idx_neutrals[n]] * m[idx_cations[c]] * pitzer.state.lambda_nc[n][c];

    // Iterate over all pairs of neutral species and anions
    for(auto n = 0; n < num_neutrals; ++n)
        for(auto a = 0; a < num_anions; ++a)
            phi += m[idx_neutrals[n]] * m[idx_anions[a]] * pitzer.state.lambda_na[n][a];

    // Iterate over all triplets of neutral species, cations and anions
    for(auto n = 0; n < num_neutrals; ++n)
        for(auto c = 0; c < num_cations; ++c)
            for(auto a = 0; a < num_anions; ++a)
                phi += m[idx_neutrals[n]] * m[idx_cations[c]] * m[idx_anions[a]] * pitzer.state.zeta[n][c][a];

    // Finalise the calculation of the osmotic coefficient
    phi = 1 + 2.0/sum_mi * phi;

    // Compute the activity of the water species
    real ln_aw = -phi * sum_mi * Mw;

    return ln_aw;
}

/// Return the Pitzer activity coefficient of a neutral species (in natural log scale).
/// @param state The state of the aqueous solution
/// @param pitzer The Pitzer parameters
/// @param N The local index of the neutral species among all neutral species in the solution
auto lnActivityCoefficientNeutral(AqueousMixtureState const& state, PitzerParams const& pitzer, Index N) -> real
{
    // The indices of the neutral species, cations and anions
    auto const& idx_cations  = pitzer.idx_cations;
    auto const& idx_anions   = pitzer.idx_anions;

    // The number of neutrals, cations and anions
    const auto num_cations  = idx_cations.size();
    const auto num_anions   = idx_anions.size();

    // The vector of molalities of all aqueous species
    auto const& m = state.m;

    // The number of species in the solution
    const auto nspecies = m.size();

    // The log of the activity coefficient of the N-th neutral species
    real ln_gammaN = {};

    // Iterate over all cations
    for(auto c = 0; c < num_cations; ++c)
        ln_gammaN += 2.0 * m[idx_cations[c]] * pitzer.state.lambda_nc[N][c];

    // Iterate over all anions
    for(auto a = 0; a < num_anions; ++a)
        ln_gammaN += 2.0 * m[idx_anions[a]] * pitzer.state.lambda_na[N][a];

    // Iterate over all pairs of cations and anions
    for(auto c = 0; c < num_cations; ++c)
        for(auto a = 0; a < num_anions; ++a)
            ln_gammaN += m[idx_cations[c]] * m[idx_anions[a]] * pitzer.state.zeta[N][c][a];

    return ln_gammaN;
}

} // namespace anonymous

auto activityModelPitzerPhreeqc(SpeciesList const& species) -> ActivityModel
{
    // Create the aqueous solution
    AqueousMixture solution(species);

    // The index of water in the solution
    const Index iwater = solution.indexWater();

    // Initialize the Pitzer params
    PitzerParams pitzer(solution);

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto stateptr = std::make_shared<AqueousMixtureState>();
    auto solutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        // Evaluate the state of the aqueous solution
        auto const& state = *stateptr = solution.state(T, P, x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // Export the aqueous solution and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = stateptr;
        props.extra["AqueousMixture"] = solutionptr;

        // Calculate the activity coefficients of the cations
        for(auto M = 0; M < pitzer.idx_cations.size(); ++M)
        {
            // Get the global index of the cation in the solution
            const Index i = pitzer.idx_cations[M];

            // Set the activity coefficient of the i-th species
            props.ln_g[i] = lnActivityCoefficientCation(state, pitzer, M);
        }

        // Calculate the activity coefficients of the anions
        for(auto X = 0; X < pitzer.idx_anions.size(); ++X)
        {
            // Get the global index of the anion in the solution
            const Index i = pitzer.idx_anions[X];

            // Set the activity coefficient of the i-th species
            props.ln_g[i] = lnActivityCoefficientAnion(state, pitzer, X);
        }

        // Calculate the activity coefficients of the neutral species
        for(auto N = 0; N < pitzer.idx_neutrals.size(); ++N)
        {
            // Get the global index of the neutral species in the solution
            const Index i = pitzer.idx_neutrals[N];

            // Set the activity coefficient of the i-th species
            props.ln_g[i] = lnActivityCoefficientNeutral(state, pitzer, N);
        }

        // Calculate the activity of water
        const real ln_aw = lnActivityWater(state, pitzer, iwater);

        // The mole fraction of water
        const auto xw = x[iwater];

        // Set the activities of the solutes
        props.ln_a = props.ln_g + log(state.m);

        // Set the activitiy of water
        props.ln_a[iwater] = ln_aw;

        // Set the activity coefficient of water (mole fraction scale)
        props.ln_g[iwater] = ln_aw - log(xw);
    };

    return fn;
}

auto ActivityModelPitzerPhreeqc() -> ActivityModelGenerator
{
    return [](SpeciesList const& species) { return activityModelPitzerPhreeqc(species); };
}

auto ActivityModelPitzer(ActivityModelParamsPitzer const& params) -> ActivityModelGenerator
{
    return [](SpeciesList const& species) { return activityModelPitzerPhreeqc(species); }; // TODO: Use params properly here!
}

auto ActivityModelPitzer(Params const& params) -> ActivityModelGenerator
{
    return [](SpeciesList const& species) { return activityModelPitzerPhreeqc(species); }; // TODO: Use params properly here!
}

} // namespace Reaktoro






























































// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright © 2014-2022 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "ActivityModelPitzer.hpp"

// // C++ includes
// // #include <set>
// // #include <string>
// // #include <vector>

// // Reaktoro includes
// #include <Reaktoro/Common/Types.hpp>
// #include <Reaktoro/Core/SpeciesList.hpp>
// // #include <Reaktoro/Common/Algorithms.hpp>
// // #include <Reaktoro/Common/ConvertUtils.hpp>
// // #include <Reaktoro/Common/Exception.hpp>
// // #include <Reaktoro/Common/NamingUtils.hpp>
// // #include <Reaktoro/Common/Real.hpp>
// // #include <Reaktoro/Common/StringUtils.hpp>
// // #include <Reaktoro/Math/BilinearInterpolator.hpp>
// // #include <Reaktoro/Water/WaterConstants.hpp>

// namespace Reaktoro {

// using std::log;
// using std::sqrt;
// using std::exp;
// using std::abs;

// enum pitz_param_type
// {
//     TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA,
//     TYPE_PSI, TYPE_ETHETA, TYPE_ALPHAS, TYPE_MU, TYPE_ETA, TYPE_Other,
//     TYPE_SIT_EPSILON, TYPE_SIT_EPSILON_MU, TYPE_APHI
// };

// struct theta_param
// {
//     real zj;
//     real zk;
//     real etheta;
//     real ethetap;
// };

// struct pitz_param
// {
//     pitz_param()
//     {
//         for(size_t i = 0; i < 3; i++) species[i] = NULL;
//         for (size_t i = 0; i < 3; i++) ispec[i] = -1;
//         type = TYPE_Other;
//         p = 0;
//         U.b0 = 0;
//         for (size_t i = 0; i < 6; i++) a[i] = 0;
//         alpha = 0;
//         os_coef = 0;
//         for (size_t i = 0; i < 3; i++) ln_coef[i] = 0;
//         thetas = NULL;
//     }
//     const char* species[3];
//     int ispec[3];
//     pitz_param_type type;
//     real p;
//     union
//     {
//         real b0;
//         real b1;
//         real b2;
//         real c0;
//         real theta;
//         real lamda;
//         real zeta;
//         real psi;
//         real alphas;
//         real mu;
//         real eta;
//         real eps;
//         real eps1;
//         real aphi;
//     } U;
//     real a[6];
//     real alpha;
//     real os_coef;
//     real ln_coef[3];
//     class theta_param* thetas;
// };

// struct PitzerPHREEQC
// {
//     int use_etheta;
//     Vec<Ptr<pitz_param>> pitz_params;
//     Vec<Ptr<theta_param>> theta_params;
//     real BK[23], DK[23];
//     real AW, VP, DW0;
// };

// auto createActivityModelPitzer(SpeciesList const& species) -> ActivityModel
// {
//     // // Inject the Pitzer namespace here
//     // using namespace Pitzer;

//     // // Create the aqueous solution
//     // AqueousMixture solution(species);

//     // // The index of water in the solution
//     // const Index iwater = solution.indexWater();

//     // // Initialize the Pitzer params
//     // PitzerParams pitzer(solution);

//     // // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
//     // auto stateptr = std::make_shared<AqueousMixtureState>();
//     // auto solutionptr = std::make_shared<AqueousMixture>(solution);

//     ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
//     {
//         // The arguments for the activity model evaluation
//         auto const& [T, P, x] = args;

//         // // Evaluate the state of the aqueous solution
//         // auto const& state = *stateptr = solution.state(T, P, x);

//         // // Set the state of matter of the phase
//         // props.som = StateOfMatter::Liquid;

//         // // Export the aqueous solution and its state via the `extra` data member
//         // props.extra["AqueousMixtureState"] = stateptr;
//         // props.extra["AqueousMixture"] = solutionptr;

//         // // Calculate the activity coefficients of the cations
//         // for(auto M = 0; M < pitzer.idx_cations.size(); ++M)
//         // {
//         //     // Get the global index of the cation in the solution
//         //     const Index i = pitzer.idx_cations[M];

//         //     // Set the activity coefficient of the i-th species
//         //     props.ln_g[i] = lnActivityCoefficientCation(state, pitzer, M);
//         // }

//         // // Calculate the activity coefficients of the anions
//         // for(auto X = 0; X < pitzer.idx_anions.size(); ++X)
//         // {
//         //     // Get the global index of the anion in the solution
//         //     const Index i = pitzer.idx_anions[X];

//         //     // Set the activity coefficient of the i-th species
//         //     props.ln_g[i] = lnActivityCoefficientAnion(state, pitzer, X);
//         // }

//         // // Calculate the activity coefficients of the neutral species
//         // for(auto N = 0; N < pitzer.idx_neutrals.size(); ++N)
//         // {
//         //     // Get the global index of the neutral species in the solution
//         //     const Index i = pitzer.idx_neutrals[N];

//         //     // Set the activity coefficient of the i-th species
//         //     props.ln_g[i] = lnActivityCoefficientNeutral(state, pitzer, N);
//         // }

//         // // Calculate the activity of water
//         // const real ln_aw = lnActivityWater(state, pitzer, iwater);

//         // // The mole fraction of water
//         // const auto xw = x[iwater];

//         // // Set the activities of the solutes
//         // props.ln_a = props.ln_g + log(state.m);

//         // // Set the activitiy of water
//         // props.ln_a[iwater] = ln_aw;

//         // // Set the activity coefficient of water (mole fraction scale)
//         // props.ln_g[iwater] = ln_aw - log(xw);
//     };

//     return fn;
// }

// auto ActivityModelPitzer() -> ActivityModelGenerator
// {
//     return [](SpeciesList const& species) { return createActivityModelPitzer(species); };
// }

// } // namespace Reaktoro
