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

/// The attributes of a substance in a fluid phase required by a cubic equation of state.
struct Substance
{
    String formula;    ///< The chemical formula of the substance in the fluid phase.
    Param Tcr = NaN;   ///< The critical temperature of the substance in the fluid phase (in K).
    Param Pcr = NaN;   ///< The critical pressure of the substance in the fluid phase (in Pa).
    Param omega = NaN; ///< The acentric factor of the substance in the fluid phase.
};

/// The properties calculated by evaluating a cubic equation of state.
struct Props
{
    real V;            ///< The molar volume of the phase (in m3/mol).
    real VT;           ///< The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).
    real VP;           ///< The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).
    real Gres;         ///< The residual molar Gibbs energy of the phase (in J/mol).
    real Hres;         ///< The residual molar enthalpy of the phase (in J/mol).
    real Cpres;        ///< The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    real Cvres;        ///< The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    ArrayXr ln_phi;    ///< The ln fugacity coefficients of the species in the phase.
    StateOfMatter som; ///< The state of matter of the fluid phase
};

/// The binary interaction parameters for the cubic equation of state.
struct Bip
{
    MatrixXr k;   ///< The computed binary interaction parameters \eq{k_{ij}} for the species in the fluid phase.
    MatrixXr kT;  ///< The computed first-order temperature derivative of the binary interaction parameters \eq{k_{ij}}.
    MatrixXr kTT; ///< The computed second-order temperature derivative of the binary interaction parameters \eq{k_{ij}}.
};

// Note: The list of arguments below for the computation of @eq{k_{ij}} and its
// first and second order temperature derivatives is inspired by the binary
// interaction parameter model of Jaubert et al. (2004):
// \eqc{k_{ij}=\dfrac{-\dfrac{1}{2}\sum_{k=1}^{N_{g}}\sum_{l=1}^{N_{g}}(\alpha_{ik}-\alpha_{jk})(\alpha_{il}-\alpha_{jl})A_{kl}\cdot(298.15T^{-1})^{\left(\frac{B_{kl}}{A_{kl}}-1\right)}-\left(\dfrac{\sqrt{a_{i}}}{b_{i}}-\dfrac{\sqrt{a_{j}}}{b_{j}}\right)}{2\dfrac{\sqrt{a_{i}a_{j}}}{b_{i}b_{j}}}.}

/// The arguments for the evaluation of binary interaction parameters in a cubic equation of state.
struct BipModelArgs
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
using BipModel = Model<Bip(BipModelArgs const&)>;

/// The value and temperature derivatives of the \eq{\alpha(T_r;\omega)} function in a cubic equation of state.
struct Alpha
{
    real alpha;   ///< The evaluated value of the \eq{\alpha(T_r;\omega)} function in a cubic equation of state.
    real alphaT;  ///< The evaluated first-order temperature derivative of the \eq{\alpha(T_r;\omega)} function.
    real alphaTT; ///< The evaluated second-order temperature derivative of the \eq{\alpha(T_r;\omega)} function.
};

/// The arguments for the evaluation of a \eq{\alpha(T_r;\omega)} function in a cubic equation of state.
struct AlphaModelArgs
{
    real const& Tr;    ///< The reduced temperature \eq{T_r=T/T_\mathrm{cr}} with respect to the substance for which \eq{\alpha(T_r;\omega)} is computed.
    real const& TrT;   ///< The first-order temperature derivative of \eq{T_r} for convenience.
    real const& omega; ///< The acentric factor of the substance for which \eq{\alpha(T_r;\omega)} is computed.
};

/// The signature of functions that evaluates \eq{\alpha(T_r;\omega)} functions for a cubic equation of state.
using AlphaModel = Model<Alpha(AlphaModelArgs const&)>;

/// The necessary constants \eq{\epsilon}, \eq{\sigma}, \eq{\Omega}, \eq{\Psi} and function \eq{\alpha(T_r;\omega)} that uniquely define a cubic equation of state.
/// We consider the following general form for a cubic equation of state \sup{\cite Smith2005}:
/// \eqc{P=\frac{RT}{V-b}-\frac{a(T)}{(V+\epsilon b)(V+\sigma b)}}
/// where:
/// \eqc{b=\sum_{i}x_{i}b_{i},}
/// \eqc{a=\sum_{i}\sum_{j}x_{i}x_{j}a_{ij},}
/// \eqc{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2},}
/// \eqc{b_{k}=\Omega\frac{RT_{\mathrm{cr},k}}{P_{\mathrm{cr},k}},}
/// \eqc{a_{k}(T)=\Psi\frac{\alpha(T_{r,k};\omega_{k})R^{2}T_{\mathrm{cr},k}^{2}}{P_{\mathrm{cr},k}}.}
/// From theq equations above, one note that a cubic equation of state can be uniquely defined
/// by constants \eq{\epsilon}, \eq{\sigma}, \eq{\Omega}, and \eq{\Psi} and function \eq{\alpha(T_r;\omega)}.
/// The table below shows how these constants and function can be defined to represent classic cubic equations of state:
///
/// | EOS | \eq{\alpha(T_{r};\omega)} | \eq{\sigma} | \eq{\epsilon} | \eq{\Omega} | \eq{\Psi} |
/// |--|--|--|--|--|--|
/// | van der Waals (1873) | 1 | 0 | 0 | 1/8 | 27/64 |
/// | Redlich-Kwong (1949)\sup{\cite Redlich1949} | \eq{T_{r}^{-1/2}} | 1 | 0 | 0.08664 | 0.42748 |
/// | Soave-Redlich-Kwong (1972) | \eq{[1+m_\mathrm{SRK}(1-\sqrt{T_{r}})]^{2}} | 1 | 0 | 0.08664 | 0.42748 |
/// | Peng-Robinson (1976)\sup{\cite Peng1976} | \eq{[1+m_\mathrm{PR76}(1-\sqrt{T_{r}})]^{2}} | \eq{1+\sqrt{2}} | \eq{1-\sqrt{2}} | 0.07780 | 0.45724 |
/// | Peng-Robinson (1978)\sup{\cite Jaubert2005} | \eq{[1+m_\mathrm{PR78}(1-\sqrt{T_{r}})]^{2}} | \eq{1+\sqrt{2}} | \eq{1-\sqrt{2}} | 0.07780 | 0.45724 |
///
/// where
/// \f[\begin{align*}m_{\mathrm{SRK}} & =0.480+1.574\omega-0.176\omega^{2}\\m_{\mathrm{PR76}} & =0.37464+1.54226\omega-0.26992\omega^{2}\vphantom{\begin{cases}\omega^2\\\omega^3\end{cases}}\\m_{\mathrm{PR78}} & =\begin{cases}0.37464+1.54226\omega-0.26992\omega^{2} & \omega\leq0.491\\0.379642+1.48503\omega-0.164423\omega^{2}+0.016666\omega^{3} & \omega>0.491\end{cases}\end{align*}\f]
struct EquationModel
{
    Param epsilon;      ///< The constant \eq{\epsilon} in the cubic equation of state.
    Param sigma;        ///< The constant \eq{\sigma} in the cubic equation of state.
    Param Omega;        ///< The constant \eq{\Omega} in the cubic equation of state.
    Param Psi;          ///< The constant \eq{\Psi} in the cubic equation of state.
    AlphaModel alphafn; ///< The function \eq{\alpha(T_r;\omega)} in the cubic equation of state.
};

/// Return a cubic equation model representative of the van der Waals (1873) cubic equation of state.
auto EquationModelVanDerWaals() -> EquationModel;

/// Return a cubic equation model representative of the Redlich-Kwong (1949) cubic equation of state.
auto EquationModelRedlichKwong() -> EquationModel;

/// Return a cubic equation model representative of the Soave-Redlich-Kwong (1972) cubic equation of state.
auto EquationModelSoaveRedlichKwong() -> EquationModel;

/// Return a cubic equation model representative of the Peng-Robinson (1976) cubic equation of state.
auto EquationModelPengRobinson76() -> EquationModel;

/// Return a cubic equation model representative of the Peng-Robinson (1978) cubic equation of state.
auto EquationModelPengRobinson78() -> EquationModel;

/// Return a cubic equation model representative of the Peng-Robinson (1978) cubic equation of state.
auto EquationModelPengRobinson() -> EquationModel;

/// The specifications for configuting a cubic equation of state.
struct EquationSpecs
{
    Vec<Substance> substances; ///< The substances in the fluid phase and their attributes.
    EquationModel eqmodel;     ///< The cubic equation of state model to be used.
    BipModel bipmodel;         ///< The function that calculates the binary interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
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

    /// Return the underlying EquationSpecs object used to create this Equation object.
    auto equationSpecs() const -> EquationSpecs const&;

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
    Param kH2O_CO2 = 0.19; ///< The binary interaction parameter \eq{k_ij} for the substance pair H₂O-CO₂.
    Param kH2O_H2S = 0.19; ///< The binary interaction parameter \eq{k_ij} for the substance pair H₂O-H₂S.
    Param kH2O_CH4 = 0.49; ///< The binary interaction parameter \eq{k_ij} for the substance pair H₂O-CH₄.
    Param kH2O_N2  = 0.49; ///< The binary interaction parameter \eq{k_ij} for the substance pair H₂O-N₂.
};

/// Return a binary interaction parameter model for Peng-Robinson EOS equivalent to that used in PHREEQC.
auto BipModelPHREEQC(Strings const& substances, BipModelParamsPHREEQC const& params = {}) -> BipModel;

} // namespace CubicEOS
} // namespace Reaktoro
