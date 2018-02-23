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

#include "TestActivityUtils.hpp"

#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

// Constants
const double nH2O  = 55.508;
const double nH    = 1e-7;
const double nOH   = 1e-7;
const double nNa   = 1.0;
const double nCl   = 1.0;
const double nNaCl = 0.2;
const double nHCl  = 1e-8;
const double nNaOH = 1e-8;

const double T     = 25.0 + 273.15;
const double P     = 1e+5;
const Vector n     = {nH2O, nH, nOH, nNa, nCl, nNaCl, nHCl, nNaOH};

const Index iH2O   = 0;
const Index iH     = 1;
const Index iOH    = 2;
const Index iNa    = 3;
const Index iCl    = 4;
const Index iNaCl  = 5;
const Index iHCl   = 6;
const Index iNaOH  = 7;

auto createAqueousMixture() -> AqueousMixture
{
    AqueousMixture mixture(8);

    mixture[0].name = "H2O(l)";
    mixture[0].charge = 0;
    mixture[0].elements = {{"H", 2}, {"O", 1}};
    mixture[0].formula = "H2O";

    mixture[1].name = "H+";
    mixture[1].charge = 1;
    mixture[1].elements = {{"H", 1}};
    mixture[1].formula = "H+";

    mixture[2].name = "OH-";
    mixture[2].charge = -1;
    mixture[2].elements = {{"H", 1}, {"O", 1}};
    mixture[2].formula = "OH-";

    mixture[3].name = "Na+";
    mixture[3].charge = 1;
    mixture[3].elements = {{"Na", 1}};
    mixture[3].formula = "Na+";

    mixture[4].name = "Cl-";
    mixture[4].charge = -1;
    mixture[4].elements = {{"Cl", 1}};
    mixture[4].formula = "Cl-";

    mixture[5].name = "NaCl(aq)";
    mixture[5].charge = 0;
    mixture[5].elements = {{"Na", 1}, {"Cl", 1}};
    mixture[5].formula = "NaCl";
    mixture[5].dissociation = {{"Na+", 1}, {"Cl-", 1}};

    mixture[6].name = "HCl(aq)";
    mixture[6].charge = 0;
    mixture[6].elements = {{"H", 1}, {"Cl", 1}};
    mixture[6].formula = "HCl";
    mixture[6].dissociation = {{"H+", 1}, {"Cl-", 1}};

    mixture[7].name = "NaOH(aq)";
    mixture[7].charge = 0;
    mixture[7].elements = {{"Na", 1}, {"O", 1}, {"H", 1}};
    mixture[7].formula = "NaOH";
    mixture[7].dissociation = {{"Na+", 1}, {"OH-", 1}};

    return mixture;
}

auto createGaseousMixture() -> GaseousMixture
{
    GaseousMixture mixture(3);

    mixture[0].name = "H2O(g)";
    mixture[0].elements = {{"H", 2}, {"O", 1}};

    mixture[1].name = "CO2(g)";
    mixture[1].elements = {{"C", 1}, {"O", 2}};

    mixture[2].name = "CH4(g)";
    mixture[2].elements = {{"C", 1}, {"H", 4}};

    return mixture;
}

auto createMineralMixture() -> MineralMixture
{
    MineralMixture mixture(2);

    mixture[0].name = "Calcite";
    mixture[0].elements = {{"Ca", 1}, {"C", 1}, {"O", 3}};

    mixture[1].name = "Magnesite";
    mixture[0].elements = {{"Mg", 1}, {"C", 1}, {"O", 3}};

    return mixture;
}

auto test_numSpecies() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(8, numSpecies(mixture));
}

auto test_indexSpecies() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, mixture.indexSpecies("H2O(l)"));
    ASSERT_EQUAL(1, mixture.indexSpecies("H+"));
    ASSERT_EQUAL(2, mixture.indexSpecies("OH-"));
    ASSERT_EQUAL(3, mixture.indexSpecies("Na+"));
    ASSERT_EQUAL(4, mixture.indexSpecies("Cl-"));
    ASSERT_EQUAL(5, mixture.indexSpecies("NaCl(aq)"));
    ASSERT_EQUAL(6, mixture.indexSpecies("HCl(aq)"));
    ASSERT_EQUAL(7, mixture.indexSpecies("NaOH(aq)"));
    ASSERT_EQUAL(numSpecies(mixture), mixture.indexSpecies(""));
    ASSERT_EQUAL(numSpecies(mixture), mixture.indexSpecies("CH4(g)"));
}

auto test_speciesNames() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    std::vector<std::string> expected = {"H2O(l)", "H+", "OH-", "Na+", "Cl-", "NaCl(aq)", "HCl(aq)", "NaOH(aq)"};
    ASSERT_EQUAL(expected, speciesNames(mixture));
}

#define ASSERT_EQUAL_ARMA_DELTA(expected, actual, delta) ASSERT(arma::norm(expected - actual)/arma::norm(actual) < delta)
#define ASSERT_EQUAL_ARMA(expected, actual) ASSERT_EQUAL_ARMA_DELTA(expected, actual, 1e-15)

auto test_speciesCharges() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Vector expected = {0, 1, -1, 1, -1, 0, 0, 0};
    ASSERT_EQUAL_ARMA(expected, speciesCharges(mixture));
}

auto test_chargedSpeciesIndices() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Indices expected = {1, 2, 3, 4};
    ASSERT_EQUAL(expected, mixture.indicesChargedSpecies());
}

auto test_chargedSpeciesLocalIndex() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, mixture.indexChargedSpecies("H+"));
    ASSERT_EQUAL(1, mixture.indexChargedSpecies("OH-"));
    ASSERT_EQUAL(2, mixture.indexChargedSpecies("Na+"));
    ASSERT_EQUAL(3, mixture.indexChargedSpecies("Cl-"));
}

auto test_chargedSpeciesNames() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    std::vector<std::string> expected = {"H+", "OH-", "Na+", "Cl-"};
    ASSERT_EQUAL(expected, chargedSpeciesNames(mixture));
}

auto test_chargedSpeciesCharges() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Vector expected = {1, -1, 1, -1};
    ASSERT_EQUAL_ARMA(expected, mixture.chargesChargedSpecies());
}

auto test_neutralSpeciesIndices() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Indices expected = {0, 5, 6, 7};
    ASSERT_EQUAL(expected, mixture.indicesNeutralSpecies());
}

auto test_neutralSpeciesLocalIndex() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, neutralSpeciesLocalIndex(mixture, "H2O(l)"));
    ASSERT_EQUAL(1, neutralSpeciesLocalIndex(mixture, "NaCl(aq)"));
    ASSERT_EQUAL(2, neutralSpeciesLocalIndex(mixture, "HCl(aq)"));
    ASSERT_EQUAL(3, neutralSpeciesLocalIndex(mixture, "NaOH(aq)"));
}

auto test_neutralSpeciesNames() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    std::vector<std::string> expected = {"H2O(l)", "NaCl(aq)", "HCl(aq)", "NaOH(aq)"};
    ASSERT_EQUAL(expected, mixture.namesNeutralSpecies());
}

auto test_cationIndices() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Indices expected = {1, 3};
    ASSERT_EQUAL(expected, mixture.indicesCations());
}

auto test_cationLocalIndex() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, cationLocalIndex(mixture, "H+"));
    ASSERT_EQUAL(1, cationLocalIndex(mixture, "Na+"));
}

auto test_cationNames() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    std::vector<std::string> expected = {"H+", "Na+"};
    ASSERT_EQUAL(expected, mixture.namesCations());
}

auto test_cationCharges() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Vector expected = {1, 1};
    ASSERT_EQUAL_ARMA(expected, mixture.chargesCations());
}

auto test_anionIndices() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Indices expected = {2, 4};
    ASSERT_EQUAL(expected, mixture.indicesAnions());
}

auto test_anionLocalIndex() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, anionLocalIndex(mixture, "OH-"));
    ASSERT_EQUAL(1, anionLocalIndex(mixture, "Cl-"));
}

auto test_anionNames() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    std::vector<std::string> expected = {"OH-", "Cl-"};
    ASSERT_EQUAL(expected, mixture.namesAnions());
}

auto test_anionCharges() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Vector expected = {-1, -1};
    ASSERT_EQUAL_ARMA(expected, mixture.chargesAnions());
}

auto test_waterIndex() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    ASSERT_EQUAL(0, mixture.indexSpecies("H2O(l)"));
}

auto test_dissociationMatrix() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    Matrix expected = zeros(4, 4);
    expected(1, 2) = 1; // expected("NaCl(aq)", "Na+") = 1
    expected(1, 3) = 1; // expected("NaCl(aq)", "Cl-") = 1
    expected(2, 0) = 1; // expected("HCl(aq)",  "H+")  = 1
    expected(2, 3) = 1; // expected("HCl(aq)",  "Cl-") = 1
    expected(3, 2) = 1; // expected("NaOH(aq)", "Na+") = 1
    expected(3, 1) = 1; // expected("NaOH(aq)", "OH-") = 1
    Matrix actual = dissociationMatrix(mixture);
    ASSERT_EQUAL_ARMA(expected, actual);
}

auto test_aqueousSolutionStateFunction() -> void
{
    AqueousMixture mixture = createAqueousMixture();
    AqueousMixtureStateFunction fn = aqueousSolutionStateFunction(mixture);
    AqueousMixtureState state = fn(T, P, n);

    auto molar_fractions_val = [=](const auto& n) -> Vector
    {
        return n/arma::sum(n);
    };

    auto molar_fractions_ddn = [](const auto& n) -> Matrix
    {
        const double sum = arma::sum(n);
        const Vector x = n/sum;
        const unsigned size = n.size();
        return arma::diagmat(x) * (arma::diagmat(1/n) - 1.0/sum*arma::ones(size, size));
    };

    auto molalities_val = [=](const auto& n) -> Vector
    {
        const double kgH2O = n[iH2O] * waterMolarMass;
        return n/kgH2O;
    };

    auto molalities_ddn = [=](const auto& n) -> Matrix
    {
        const double kgH2O = n[iH2O] * waterMolarMass;
        const unsigned size = n.size();
        Matrix ddn = arma::eye(size, size);
        ddn.col(iH2O) -= n/n[iH2O];
        return ddn/kgH2O;
    };

    VectorFunction stoichiometric_molalities_fn = [=](const auto& n)
    {
        const double kgH2O = n[iH2O] * waterMolarMass;
        const double msH  = 1.0/kgH2O * (n[iH]  + n[iHCl]);
        const double msOH = 1.0/kgH2O * (n[iOH] + n[iNaOH]);
        const double msNa = 1.0/kgH2O * (n[iNa] + n[iNaCl] + n[iNaOH]);
        const double msCl = 1.0/kgH2O * (n[iCl] + n[iNaCl] + n[iHCl]);
        return Vector{msH, msOH, msNa, msCl};
    };

    ScalarFunction effective_ionic_strengh_fn = [=](const auto& n)
    {
        Vector m = molalities_val(n);
        return 0.5 * (m[iH] + m[iOH] + m[iNa] + m[iCl]);
    };

    ScalarFunction stoichiometric_ionic_strengh_fn = [=](const auto& n)
    {
        Vector ms = stoichiometric_molalities_fn(n);
        return 0.5 * (ms[0] + ms[1] + ms[2] + ms[3]); // use local indices
    };

    Vector x_val  = molar_fractions_val(n);
    Vector m_val  = molalities_val(n);
    Vector ms_val = stoichiometric_molalities_fn(n);
    double Ie_val = effective_ionic_strengh_fn(n);
    double Is_val = stoichiometric_ionic_strengh_fn(n);

    Matrix x_ddn  = molar_fractions_ddn(n);
    Matrix m_ddn  = molalities_ddn(n);
    Matrix ms_ddn = derivativeForward(stoichiometric_molalities_fn, n);
    Vector Ie_ddn = derivativeForward(effective_ionic_strengh_fn, n);
    Vector Is_ddn = derivativeForward(stoichiometric_ionic_strengh_fn, n);

    ASSERT_EQUAL(T, state.T);
    ASSERT_EQUAL(P, state.P);
    ASSERT_EQUAL_ARMA(n, state.n);
    ASSERT_EQUAL_DELTA(Ie_val, state.Ie.val, 1e-15);
    ASSERT_EQUAL_DELTA(Is_val, state.Is.val, 1e-15);
    ASSERT_EQUAL_ARMA(x_val, state.x.val);
    ASSERT_EQUAL_ARMA(m_val, state.m.val);
    ASSERT_EQUAL_ARMA(ms_val, state.ms.val);

    ASSERT_EQUAL_ARMA_DELTA(x_ddn, state.x.ddn, 1e-6);
    ASSERT_EQUAL_ARMA_DELTA(m_ddn, state.m.ddn, 1e-6);
    ASSERT_EQUAL_ARMA_DELTA(ms_ddn, state.ms.ddn, 1e-6);
    ASSERT_EQUAL_ARMA_DELTA(Ie_ddn, state.Ie.ddn, 1e-6);
    ASSERT_EQUAL_ARMA_DELTA(Is_ddn, state.Is.ddn, 1e-6);
}

auto test_gaseousSolutionStateFunction() -> void
{
    GaseousMixture mixture = createGaseousMixture();
    GaseousMixtureStateFunction fn = gaseousSolutionStateFunction(mixture);
    GaseousMixtureState state = fn(T, P, n);

    VectorFunction molar_fractions_fn = [=](const auto& n)
    {
        return n/arma::sum(n);
    };

    Vector x_val  = molar_fractions_fn(n);
    Matrix x_ddn  = derivativeCentral(molar_fractions_fn, n);

    ASSERT_EQUAL(T, state.T);
    ASSERT_EQUAL(P, state.P);
    ASSERT_EQUAL_ARMA(n, state.n);
    ASSERT_EQUAL_ARMA(x_val, state.x.val);
}

auto test_mineralSolutionStateFunction() -> void
{
    MineralMixture mixture = createMineralMixture();
    MineralMixtureStateFunction fn = mineralSolutionStateFunction(mixture);
    MineralMixtureState state = fn(T, P, n);

    VectorFunction molar_fractions_fn = [=](const auto& n)
    {
        return n/arma::sum(n);
    };

    Vector x_val  = molar_fractions_fn(n);
    Matrix x_ddn  = derivativeCentral(molar_fractions_fn, n);

    ASSERT_EQUAL(T, state.T);
    ASSERT_EQUAL(P, state.P);
    ASSERT_EQUAL_ARMA(n, state.n);
    ASSERT_EQUAL_ARMA(x_val, state.x.val);
}

} // namespace

auto testSuiteActivityUtils() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_numSpecies);
    s += CUTE(test_indexSpecies);
    s += CUTE(test_speciesNames);
    s += CUTE(test_speciesCharges);
    s += CUTE(test_chargedSpeciesIndices);
    s += CUTE(test_chargedSpeciesLocalIndex);
    s += CUTE(test_chargedSpeciesNames);
    s += CUTE(test_chargedSpeciesCharges);
    s += CUTE(test_neutralSpeciesIndices);
    s += CUTE(test_neutralSpeciesLocalIndex);
    s += CUTE(test_neutralSpeciesNames);
    s += CUTE(test_cationIndices);
    s += CUTE(test_cationLocalIndex);
    s += CUTE(test_cationNames);
    s += CUTE(test_cationCharges);
    s += CUTE(test_anionIndices);
    s += CUTE(test_anionLocalIndex);
    s += CUTE(test_anionNames);
    s += CUTE(test_anionCharges);
    s += CUTE(test_waterIndex);
    s += CUTE(test_dissociationMatrix);
    s += CUTE(test_aqueousSolutionStateFunction);
    s += CUTE(test_gaseousSolutionStateFunction);
    s += CUTE(test_mineralSolutionStateFunction);

    return s;
}

} // namespace Reaktoro
