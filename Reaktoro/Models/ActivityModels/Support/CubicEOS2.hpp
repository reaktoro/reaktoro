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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {
namespace CubicEOS {

/// The properties calculated by evaluating a cubic equation of state.
struct Props
{
    /// The molar volume of the phase (in m3/mol).
    real V;

    /// The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).
    real VT;

    /// The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).
    real VP;

    /// The residual molar Gibbs energy of the phase (in J/mol).
    real Gres;

    /// The residual molar enthalpy of the phase (in J/mol).
    real Hres;

    /// The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    real Cpres;

    /// The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    real Cvres;

    /// The ln fugacity coefficients of the species in the phase.
    ArrayXr ln_phi;

    /// The state of matter of the fluid phase
    StateOfMatter som;
};

/// The currently supported types of cubic equations of state.
enum class Type
{
    VanDerWaals,
    RedlichKwong,
    SoaveRedlichKwong,
    PengRobinson
};

/// The binary interaction parameters for the cubic equation of state.
struct Bip
{
    /// The computed binary interaction parameters \eq{k_{ij}} for the species in the fluid phase.
    MatrixXr k;

    /// The computed first-order temperature derivative of the binary interaction parameters \eq{k_{T,ij}}.
    MatrixXr kT;

    /// The computed second-order temperature derivative of the binary interaction parameters \eq{k_{TT,ij}}.
    MatrixXr kTT;
};

// Note: The list of arguments below for the computation of @eq{k_{ij}} and its
// first and second order temperature derivatives is inspired by the binary
// interaction parameter model of Jaubert et al. (2004):
// \eqc{k_{ij}=\dfrac{-\dfrac{1}{2}\sum_{k=1}^{N_{g}}\sum_{l=1}^{N_{g}}(\alpha_{ik}-\alpha_{jk})(\alpha_{il}-\alpha_{jl})A_{kl}\cdot(298.15T^{-1})^{\left(\frac{B_{kl}}{A_{kl}}-1\right)}-\left(\dfrac{\sqrt{a_{i}}}{b_{i}}-\dfrac{\sqrt{a_{j}}}{b_{j}}\right)}{2\dfrac{\sqrt{a_{i}a_{j}}}{b_{i}b_{j}}}.}

/// The critical properties of a substance in the fluid phase.
struct CrProps
{
    Param Tcr   = NaN; ///< The critical temperature of the substance in the fluid phase (in K).
    Param Pcr   = NaN; ///< The critical pressure of the substance in the fluid phase (in Pa).
    Param omega = NaN; ///< The acentric factor of the substance in the fluid phase.
};

/// The arguments for the evaluation of binary interaction parameters in a cubic equation of state.
struct BipArgs
{
    Strings const& substances; ///< The chemical formulas of the substances in the fluid phase.
    real const& T;             ///< The temperature of the fluid (in K).
    ArrayXrConstRef Tcr;       ///< The critical temperatures of the substances in the fluid phase (in K).
    ArrayXrConstRef Pcr;       ///< The critical pressures of the substances in the fluid phase (in Pa).
    ArrayXrConstRef omega;     ///< The acentric factors of the substances in the fluid phase.
    ArrayXrConstRef a;         ///< The variables \eq{a_i} in the cubic equation of state for each substance in the fluid phase computed at current temperature.
    ArrayXrConstRef aT;        ///< The first-order temperature derivative of \eq{a_i}.
    ArrayXrConstRef aTT;       ///< The second-order temperature derivative of \eq{a_i}.
    ArrayXrConstRef alpha;     ///< The variables \eq{\alpha_i} in the cubic equation of state for each substance in the fluid phase computed at current temperature.
    ArrayXrConstRef alphaT;    ///< The first-order temperature derivative of \eq{\alpha_i}.
    ArrayXrConstRef alphaTT;   ///< The second-order temperature derivative of \eq{\alpha_i}.
    ArrayXrConstRef b;         ///< The variables \eq{b_i} in the cubic equation of state for each substance in the fluid phase computed at current temperature.
};

/// The signature of functions that evaluates interaction parameters for a cubic equation of state.
using BipModel = Model<Bip(BipArgs const&)>;

struct Alpha
{
    real alpha;
    real alphaT;
    real alphaTT;
};

struct AlphaModelArgs
{
    real const& Tr;
    real const& TrT;
    real const& omega;
};

using AlphaModel = Model<Alpha(AlphaModelArgs const&)>;

struct EquationModel
{
    Param sigma;
    Param epsilon;
    Param Omega;
    Param Psi;
    AlphaModel alphafn;
};

auto EquationModelVanDerWaals() -> EquationModel;
auto EquationModelRedlichKwong() -> EquationModel;
auto EquationModelSoaveRedlichKwong() -> EquationModel;
auto EquationModelPengRobinson() -> EquationModel;

/// The specifications for configuting a cubic equation of state.
struct EquationSpecs
{
    /// The substances in the fluid phase and their critical properties.
    Dict<String, CrProps> substances;

    /// The cubic equation of state model to be used.
    EquationModel eqmodel = EquationModelPengRobinson();

    /// The function that calculates interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
    BipModel bipmodel;
};

/// Calculates thermodynamic properties of fluid phases based on a cubic equation of state model.
class Equation
{
public:
    /// Construct an Equation object representing a cubic equation of state.
    explicit Equation(EquationSpecs const& specs);

    /// Construct a copy of an Equation object.
    Equation(Equation const& other);

    /// Destroy this Equation object.
    virtual ~Equation();

    /// Assign an Equation object to this.
    auto operator=(Equation other) -> Equation&;

    /// Compute the thermodynamic properties of the phase.
    /// @param[in] props The evaluated thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in K)
    /// @param P The pressure of the phase (in Pa)
    /// @param x The mole fractions of the species in the phase (in mol/mol)
    auto compute(Props& props, real const& T, real const& P, ArrayXrConstRef const& x) -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// The parameters in the binary interaction parameter model of PHREEQC for Peng-Robinson EOS.
/// The parameters below were collected from method `Phreeqc::calc_PR` in file `gases.cpp` from PHREEQC source code.
struct BipModelParamsPHREEQC
{
    Param kH2O_CO2 = 0.19;
    Param kH2O_H2S = 0.19;
    Param kH2O_CH4 = 0.49;
    Param kH2O_N2  = 0.49;
};

auto BipModelPHREEQC(Strings const& substances, BipModelParamsPHREEQC const& params = {}) -> BipModel;

} // namespace CubicEOS
} // namespace Reaktoro
