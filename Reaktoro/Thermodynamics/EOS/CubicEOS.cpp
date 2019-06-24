// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "CubicEOS.hpp"

// C++ includes
#include <algorithm>
#include <iostream>  // TODO: REMOVE
#include <limits>
#include <Reaktoro/Math/Matrix.hpp> // TODO: REMOVE
#include <fstream>
#include <array>
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>

#include <complex>

#include "W:\release\Projects\Reaktoro\demos\cpp\phaseid.hpp"

#define M_PI 3.14159265358979323846

extern phaseIdentificationMethod phaseidMethod;

int quantity = 0; 

namespace Reaktoro {
namespace internal {

using AlphaResult = std::tuple<ThermoScalar, ThermoScalar, ThermoScalar>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT (temperature derivatives) for a given EOS.
auto alpha(CubicEOS::Model type) -> std::function<AlphaResult(const ThermoScalar&, double)>
{
    // The alpha function for van der Waals EOS
    auto alphaVDW = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val(1.0);
        ThermoScalar ddt(0.0);
        ThermoScalar d2dt2(0.0);
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val = 1.0/sqrt(Tr);
        ThermoScalar ddt = -0.5/Tr * val;
        ThermoScalar d2dt2 = -0.5/Tr * (ddt - val/Tr);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        double m = 0.480 + 1.574*omega - 0.176*omega*omega;
        ThermoScalar sqrtTr = sqrt(Tr);
        ThermoScalar aux_val = 1.0 + m*(1.0 - sqrtTr);
        ThermoScalar aux_ddt = -0.5*m/sqrtTr;
        ThermoScalar aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        ThermoScalar val = aux_val*aux_val;
        ThermoScalar ddt = 2.0*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.
        double m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        ThermoScalar sqrtTr = sqrt(Tr);
        ThermoScalar aux_val = 1.0 + m*(1.0 - sqrtTr);
        ThermoScalar aux_ddt = -0.5*m/sqrtTr;
        ThermoScalar aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        ThermoScalar val = aux_val*aux_val;
        ThermoScalar ddt = 2.0*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    switch(type)
    {
        case CubicEOS::VanDerWaals: return alphaVDW;
        case CubicEOS::RedlichKwong: return alphaRK;
        case CubicEOS::SoaveRedlichKwong: return alphaSRK;
        case CubicEOS::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto sigma(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 1.0;
        case CubicEOS::SoaveRedlichKwong: return 1.0;
        case CubicEOS::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
    }
}

auto epsilon(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 0.0;
        case CubicEOS::SoaveRedlichKwong: return 0.0;
        case CubicEOS::PengRobinson: return 1.0 - 1.4142135623730951;
        default: return 1.0 - 1.4142135623730951;
    }
}

auto Omega(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 1.0/8.0;
        case CubicEOS::RedlichKwong: return 0.08664;
        case CubicEOS::SoaveRedlichKwong: return 0.08664;
        case CubicEOS::PengRobinson: return 0.0777960739;
        default: return 0.0777960739;
    }
}

auto Psi(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 27.0/64.0;
        case CubicEOS::RedlichKwong: return 0.42748;
        case CubicEOS::SoaveRedlichKwong: return 0.42748;
        case CubicEOS::PengRobinson: return 0.457235529;
        default: return 0.457235529;
    }
}

} // namespace internal


namespace GibbsEnergyCriticPressure {
	auto should_be_gas() -> bool {
		return true;
	}
}

namespace VolumeMethod {
	auto should_be_gas(double V, double bmix) -> bool {
		return (V/bmix) > 1.75;
	}
}

namespace CriticalPointMethods {
	auto should_be_gas(double amix, double bmix, double R, double Omega_b, double Omega_a, double Z, double V, double T) -> bool {
		auto Tc = (1.0 / R)*(amix / bmix)*(Omega_b / Omega_a);
		auto Vc = Z * (bmix / Omega_b);
		return (V * T*T) > (Vc*Tc*Tc);
	}
}

namespace IsothermalCompressibilityMethods {
	auto should_be_gas(ThermoScalar T, ThermoScalar P, ChemicalScalar Z, double R) -> bool {
		/*
		double dT = 0.01;
		auto T_minos = T - dT;
		auto T_plus = T + dT;
		auto V_minos = Z * R*T_minos / P;
		auto V_plus = Z * R*T_plus / P;
		auto k_minos = -(1.0 / V_minos.val) * V_minos.ddP;
		auto k_plus = -(1.0 / V_plus.val) * V_plus.ddP;
		auto dkdT = (k_plus - k_minos) / dT;
		
		auto V = Z * R * T / P;
		auto k = -(1.0 / V.val) * V.ddP;
		*/

		auto V = Z * R * T / P;

		
		auto dkdt = (1.0 / (V.val*V.val))*V.ddP*V.ddT;
		


		return  dkdt <= 0.0; // k*P.val >= 0.9 && k*P.val <= 3.0;// true;//
	}
}

namespace workanalisys {
	auto should_be_gas(double Vgas, double Vliq, double a, double b, double R, double T, double P) -> bool {

		auto w1 = P * (Vgas - Vliq);
		auto w2 = R * T*std::log((Vgas - b) / (Vliq - b)) + (a / (std::sqrt(T)*b))*std::log(((Vgas + b)*Vliq) / ((Vliq + b)*Vgas));

		return w2 - w1 > 0;
	};
}

namespace workanalisys_PengRobinson {
	auto should_be_gas(double Vgas, double Vliq, double a, double b, double R, double T, double P) -> bool {

		auto w1 = P * (Vgas - Vliq);
		auto w2_a = R * T*std::log((Vgas - b) / (Vliq - b));
		auto aux_gas = (b + Vgas) / (2 * std::sqrt(2)*b);
		auto aux_liq = (b + Vgas) / (2 * std::sqrt(2)*b);
		auto w2_b = a * ((std::log(1 - ((b + Vgas) / (2 * std::sqrt(2)*b))) / (2 * std::sqrt(2)*b)) - (std::log(((b + Vgas) / (2 * std::sqrt(2)*b)) + 1) / (2 * std::sqrt(2)*b)) -
			(std::log(1 - ((b + Vliq) / (2 * std::sqrt(2)*b))) / (2 * std::sqrt(2)*b)) + (std::log(((b + Vliq) / (2 * std::sqrt(2)*b)) + 1) / (2 * std::sqrt(2)*b)));
		
		
		auto w2_c = a * ((std::log(1 - aux_gas) / (std::sqrt(2)*b)) - (std::log(aux_gas + 1) / (std::sqrt(2)*b)) -
			(std::log(1 - aux_liq) / (std::sqrt(2)*b)) + (std::log(aux_liq + 1) / (std::sqrt(2)*b)));
		
		//auto w2 = R * T*std::log((Vgas - b) / (Vliq - b)) + (a / (std::sqrt(T)*b))*std::log(((Vgas + b)*Vliq) / ((Vliq + b)*Vgas));

		if ((w2_a + w2_b - w1) > 0.0) {
			std::cout << "true" << std::endl;
		}

		return w2_a + w2_b - w1 > 0;
	};
}

struct CubicEOS::Impl
{
    /// The number of species in the phase.
    unsigned nspecies;

    /// The flag that indicates if the phase is vapor (false means liquid instead).
    bool isvapor = true;

	bool minimaze_gibbs_energy = false;

    /// The type of the cubic equation of state.
    CubicEOS::Model model = CubicEOS::PengRobinson;

    /// The critical temperatures of the species (in units of K).
    std::vector<double> critical_temperatures;

    /// The critical pressures of the species (in units of Pa).
    std::vector<double> critical_pressures;

    /// The acentric factors of the species.
    std::vector<double> acentric_factors;

    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    InteractionParamsFunction calculate_interaction_params;

    /// The result with thermodynamic properties calculated from the cubic equation of state
    Result result;

    /// Construct a CubicEOS::Impl instance.
    Impl(unsigned nspecies)
    : nspecies(nspecies)
    {
        // Initialize the dimension of the chemical vector quantities
        ChemicalVector vec(nspecies);
        result.partial_molar_volumes = vec;
        result.residual_partial_molar_enthalpies = vec;
        result.residual_partial_molar_gibbs_energies = vec;
        result.ln_fugacity_coefficients = vec;
    }

	
	
		auto unique_root (std::complex<double> a, std::complex<double> b, std::complex<double> c) -> bool {
			if (a == b) {
				return false;
			}
			else if (b == c) {
				return false;
			}
			else if (a == c) {
				return false;
			}
			else {
				return true;
			}
		}

		auto there_is_img (std::complex<double> a, std::complex<double> b, std::complex<double> c) -> bool {
			if (a.imag() != 0 || b.imag() != 0 || c.imag() != 0) {
				return true;
			}
			else {
				return false;
			}
		}

		auto use_phase_identivation (CubicRoots roots) -> bool {


			return !(unique_root(std::get<0>(roots), std::get<1>(roots), std::get<2>(roots))) && there_is_img(std::get<0>(roots), std::get<1>(roots), std::get<2>(roots));
		}

	
		auto get_real_roots(std::complex<double> a, std::complex<double> b, std::complex<double> c) -> std::vector<double>{
			std::vector<double> real_roots;
			if (a.imag() == 0) {
				real_roots.push_back(a.real());
			}
			if (b.imag() == 0) {
				real_roots.push_back(b.real());
			}
			if (c.imag() == 0) {
				real_roots.push_back(c.real());
			}
			return real_roots;
		}
	



    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> Result
    {
		quantity++;
        // Check if the mole fractions are zero or non-initialized
        if(x.val.size() == 0 || min(x.val) <= 0.0)
            return Result(nspecies); // result with zero values

        // Auxiliary variables
        const double R = universalGasConstant;
        const double Psi = internal::Psi(model);
        const double Omega = internal::Omega(model);
        const double epsilon = internal::epsilon(model);
        const double sigma = internal::sigma(model);
        const auto alpha = internal::alpha(model);

        // Calculate the parameters `a` of the cubic equation of state for each species
        ThermoVector a(nspecies);
        ThermoVector aT(nspecies);
        ThermoVector aTT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            const double factor = Psi*R*R*(Tc*Tc)/Pc;
            const ThermoScalar Tr = T/Tc;
            ThermoScalar alpha_val, alpha_ddt, alpha_d2dt2;
            std::tie(alpha_val, alpha_ddt, alpha_d2dt2) = alpha(Tr, omega);
            a[i] = factor * alpha_val;
            aT[i] = factor * alpha_ddt;
            aTT[i] = factor * alpha_d2dt2;
        };

        // Calculate the parameters `b` of the cubic equation of state for each species
        Vector b(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            b[i] = Omega*R*Tci/Pci;
        }

        // Calculate the table of binary interaction parameters and its temperature derivatives
        InteractionParamsResult kres;
        InteractionParamsArgs kargs{T, a, aT, aTT, b};

        if(calculate_interaction_params)
            kres = calculate_interaction_params(kargs);
		
		/*
		std::vector<ThermoScalar> H2S = { ThermoScalar(0.0000), ThermoScalar(0.0582) };
		std::vector<ThermoScalar> CH4 = { ThermoScalar(0.0582), ThermoScalar(0.0000) };
		*/

		/*
		std::vector<ThermoScalar> CH4_g = { ThermoScalar(0.0000)								, ThermoScalar( 0.0497), ThermoScalar(0.0582), ThermoScalar(-5.33e-6*T.val*T.val + 4.73e-3*T.val - 0.947) };
		std::vector<ThermoScalar> CO2_g = { ThermoScalar(0.0497)								, ThermoScalar( 0.0000), ThermoScalar(0.0669), ThermoScalar(-0.0169)									  };
		std::vector<ThermoScalar> H2S_g = { ThermoScalar(0.0582)								, ThermoScalar( 0.0669), ThermoScalar(0.0000), ThermoScalar( 0.0362)									  };
		std::vector<ThermoScalar> H2O_g = { ThermoScalar(-5.33e-6*T.val + 4.73e-3*T.val - 0.947), ThermoScalar(-0.0169), ThermoScalar(0.0362), ThermoScalar( 0.0000)									  };
		
		std::vector<ThermoScalar> CH4_gT = { ThermoScalar(0.0000)					, ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(-10.66e-6*T.val + 4.73e-3) };
		std::vector<ThermoScalar> CO2_gT = { ThermoScalar(0.0000)				    , ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)					  };
		std::vector<ThermoScalar> H2S_gT = { ThermoScalar(0.0000)				    , ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)					  };
		std::vector<ThermoScalar> H2O_gT = { ThermoScalar(-10.66e-6*T.val + 4.73e-3), ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)					  };

		std::vector<ThermoScalar> CH4_gTT = { ThermoScalar(0.0000)	 , ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(-10.66e-6) };
		std::vector<ThermoScalar> CO2_gTT = { ThermoScalar(0.0000)	 , ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)	   };
		std::vector<ThermoScalar> H2S_gTT = { ThermoScalar(0.0000)	 , ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)    };
		std::vector<ThermoScalar> H2O_gTT = { ThermoScalar(-10.66e-6), ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000)    };
		*/
		/*
		
		std::vector<ThermoScalar> CH4_oil = { ThermoScalar(0.0000), ThermoScalar(0.0497), ThermoScalar(0.0582) };
		std::vector<ThermoScalar> CO2_oil = { ThermoScalar(0.0497), ThermoScalar(0.0000), ThermoScalar(0.0669) };
		std::vector<ThermoScalar> H2S_oil = { ThermoScalar(0.0582), ThermoScalar(0.0669), ThermoScalar(0.0000) };

		std::vector<ThermoScalar> CH4_oilT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };
		std::vector<ThermoScalar> CO2_oilT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };
		std::vector<ThermoScalar> H2S_oilT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };

		std::vector<ThermoScalar> CH4_oilTT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };
		std::vector<ThermoScalar> CO2_oilTT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };
		std::vector<ThermoScalar> H2S_oilTT = { ThermoScalar(0.0000), ThermoScalar(0.0000), ThermoScalar(0.0000) };
		*/
		/*
		//if (isvapor) {
			kres.k.push_back(CH4_g);
			kres.k.push_back(CO2_g);
			kres.k.push_back(H2S_g);
			kres.k.push_back(H2O_g);
			kres.kT.push_back(CH4_gT);
			kres.kT.push_back(CO2_gT);
			kres.kT.push_back(H2S_gT);
			kres.kT.push_back(H2O_gT);
			kres.kTT.push_back(CH4_gTT);
			kres.kTT.push_back(CO2_gTT);
			kres.kTT.push_back(H2S_gTT);
			kres.kTT.push_back(H2O_gTT);
		//}
		*/

		
			/*
		if (!isvapor) {
			kres.k.push_back(CH4_oil);
			kres.k.push_back(CO2_oil);
			kres.k.push_back(H2S_oil);
			kres.kT.push_back(CH4_oilT);
			kres.kT.push_back(CO2_oilT);
			kres.kT.push_back(H2S_oilT);
			kres.kTT.push_back(CH4_oilTT);
			kres.kTT.push_back(CO2_oilTT);
			kres.kTT.push_back(H2S_oilTT);

		}
		*/
		//kres.k.push_back(H2S);
		//kres.k.push_back(CH4);		

        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        ChemicalScalar amix(nspecies);
        ChemicalScalar amixT(nspecies);
        ChemicalScalar amixTT(nspecies);
        ChemicalVector abar(nspecies);
        ChemicalVector abarT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const ThermoScalar r = kres.k.empty() ? ThermoScalar(1.0) : 1.0 - kres.k[i][j];
                const ThermoScalar rT = kres.kT.empty() ? ThermoScalar(0.0) : -kres.kT[i][j];
                const ThermoScalar rTT = kres.kTT.empty() ? ThermoScalar(0.0) : -kres.kTT[i][j];

                const ThermoScalar s = sqrt(a[i]*a[j]);
                const ThermoScalar sT = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const ThermoScalar sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const ThermoScalar aij = r*s;
                const ThermoScalar aijT = rT*s + r*sT;
                const ThermoScalar aijTT = rTT*s + 2.0*rT*sT + r*sTT;

                amix += x[i] * x[j] * aij;
                amixT += x[i] * x[j] * aijT;
                amixTT += x[i] * x[j] * aijTT;

                abar[i] += 2 * x[j] * aij;
                abarT[i] += 2 * x[j] * aijT;
            }
        }

        // Finalize the calculation of `abar` and `abarT`
        for(unsigned i = 0; i < nspecies; ++i)
        {
            abar[i] -= amix;
            abarT[i] -= amixT;
        }

        // Calculate the parameter `bmix` of the cubic equation of state
        ChemicalScalar bmix(nspecies);
        Vector bbar(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            bbar[i] = Omega*R*Tci/Pci;
            bmix += x[i] * bbar[i];
        }

        // Calculate the temperature derivative of `bmix`
        const double bmixT = 0.0; // no temperature dependence

        // Calculate auxiliary quantities `beta` and `q`
        const ChemicalScalar beta = P*bmix/(R*T);
        const ChemicalScalar betaT = beta * (bmixT/bmix - 1.0/T);

        const ChemicalScalar q = amix/(bmix*R*T);
        const ChemicalScalar qT = q*(amixT/amix - 1.0/T);
        const ChemicalScalar qTT = qT*qT/q + q*(1.0/(T*T) + amixTT/amix - amixT*amixT/(amix*amix));

        // Calculate the coefficients A, B, C of the cubic equation of state
        const ChemicalScalar A = (epsilon + sigma - 1)*beta - 1;
        const ChemicalScalar B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
        const ChemicalScalar C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate the partial temperature derivative of the coefficients A, B, C
        const ChemicalScalar AT = (epsilon + sigma - 1)*betaT;
        const ChemicalScalar BT = 2*(epsilon*sigma - epsilon - sigma)*beta*betaT + qT*beta - (epsilon + sigma - q)*betaT;
        const ChemicalScalar CT = -3*epsilon*sigma*beta*beta*betaT - qT*beta*beta - 2*(epsilon*sigma + q)*beta*betaT;

        // Define the non-linear function and its derivative for calculation of its root
        const auto f = [&](double Z) -> std::tuple<double, double>
        {
            const double val = Z*Z*Z + A.val*Z*Z + B.val*Z + C.val;
            const double grad = 3*Z*Z + 2*A.val*Z + B.val;
            return std::make_tuple(val, grad);
        };

        // Define the parameters for Newton's method
        const auto tolerance = 1e-6;
        const auto maxiter = 10000000;

        // Determine the appropriate initial guess for the cubic equation of state
        const double Z0 = isvapor ? 1.0 : beta.val;

        // Calculate the compressibility factor Z using Newton's method
        ChemicalScalar Z(nspecies);

        //auto roots_result = calculate_roots(A.val, B.val, C.val);
        //auto has_real_roots = std::get<0>(roots_result);
        //auto smallest_root = std::get<1>(roots_result);
        //auto greatest_root = std::get<2>(roots_result);
        
        //std::cout << "cardano method" << std::endl;
        auto roots_im = cardano(1, A.val, B.val, C.val);
        auto roots = { std::get<0>(roots_im), std::get<1>(roots_im), std::get<2>(roots_im) };
        auto first_root = std::get<0>(roots_im);
        auto second_root = std::get<1>(roots_im);
        auto third_root = std::get<2>(roots_im);
        auto modulo1 = std::sqrt(first_root.real()*first_root.real() + first_root.imag()*first_root.imag());
        auto modulo2 = std::sqrt(second_root.real()*second_root.real() + second_root.imag()*second_root.imag());
        auto modulo3 = std::sqrt(third_root.real()*third_root.real() + third_root.imag()*third_root.imag());
        //std::cout << "imt_root1" << first_root << std::endl;
        //std::cout << "modulo1" << std::sqrt(first_root.real()*first_root.real() + first_root.imag()*first_root.imag()) << std::endl;
        //std::cout << "imt_root2" << second_root << std::endl;
        //std::cout << "modulo2" << std::sqrt(second_root.real()*second_root.real() + second_root.imag()*second_root.imag()) << std::endl;
        //std::cout << "imt_root3" << third_root << std::endl;
        //std::cout << "modulo3" << std::sqrt(third_root.real()*third_root.real() + third_root.imag()*third_root.imag()) << std::endl;
        //Z.val = newton(f, Z0, tolerance, maxiter);
       
		auto real_roots = get_real_roots(first_root, second_root, third_root);

		//auto teste = use_phase_identivation(roots_im);

		/*
		if (first_root.real() < 0) {
			first_root.real() = 0;
		}
		*/

        auto comp1 = [](auto a, auto b) {
            return a.real() > b.real();
        };

        auto comp2 = [](auto a, auto b) {
            return std::sqrt(a.real()*a.real() + a.imag()*a.imag()) > std::sqrt(b.real()*b.real() + b.imag()*b.imag());
        };

		auto comp3 = [](auto a, auto b) {
			return a > b;
		};

		if (isvapor) {
			Z.val = *std::max_element(std::begin(real_roots), std::end(real_roots));
		}
		else {
			Z.val = *std::min_element(std::begin(real_roots), std::end(real_roots));
		}

		ChemicalScalar& V = result.molar_volume;
		ChemicalScalar& G_res = result.residual_molar_gibbs_energy;
		ChemicalScalar& H_res = result.residual_molar_enthalpy;
		ChemicalScalar& Cp_res = result.residual_molar_heat_capacity_cp;
		ChemicalScalar& Cv_res = result.residual_molar_heat_capacity_cv;
		ChemicalVector& Vi = result.partial_molar_volumes;
		ChemicalVector& Gi_res = result.residual_partial_molar_gibbs_energies;
		ChemicalVector& Hi_res = result.residual_partial_molar_enthalpies;
		ChemicalVector& ln_phi = result.ln_fugacity_coefficients;

        
		bool should_be_gas = false;

		bool should_be_gas_gibbs = false;

		auto Zs_pair = std::minmax_element(std::begin(real_roots), std::end(real_roots));

		std::array<ChemicalScalar, 2> Zs;
		Zs[0] = ChemicalScalar(nspecies, *(Zs_pair).first);
		Zs[1] = ChemicalScalar(nspecies, *(Zs_pair).second);
		auto V_liq = R * Zs[0] * T / P;
		auto V_gas = R * Zs[1] * T / P;

		if (real_roots.size() == 1) {
			const auto p = [&](double V) -> double
			{
				const double val = ((R*T.val) / (V - bmix.val)) - (amix.val / ((V + epsilon * bmix.val) * (V + sigma * bmix.val)));
				return val;
			};
			const auto dp_dv = [&](double V) -> double
			{
				double k1 = V + epsilon * bmix.val;
				double k2 = V + sigma * bmix.val;
				const double val = -(R * T.val) / (std::pow(V - bmix.val, 2.0)) + amix.val * (k1 + k2) / (std::pow(k1, 2.0) * std::pow(k2, 2.0));
				return val;
			};
			auto real_cubic_root = [](double x) -> double
			{
				if (x < 0.0) {
					return -(std::pow(std::abs(x), 1.0 / 3.0));
				}
				else if (x == 0) {
					return 0.0;
				}
				else {
					return (std::pow(std::abs(x), 1.0 / 3.0));
				}
			};

			auto cubic_roots = [real_cubic_root](double a, double b, double c, double d) -> std::tuple<std::complex<double>, std::complex<double>, std::complex<double>>
			{
				//a x ^3 + b x^2 + c x + d =0
				auto p = c / a - 1.0 / 3.0 * std::pow(b / a, 2.0);
				auto q = (2.0 / 27.0) * std::pow(b / a, 3.0) - b * c / (3.0 * std::pow(a, 2.0)) + d / a;
				auto delta = std::pow(q, 2.0) / 4.0 + std::pow(p, 3.0) / 27.0;

				if (delta > 0) {
					auto u = real_cubic_root(-q / 2.0 + std::pow(delta, 1.0 / 2.0));
					auto v = real_cubic_root(-q / 2.0 - std::pow(delta, 1.0 / 2.0));
					std::complex<double> x1, x2, x3;
					x1.real(-b / (3.0 * a) + u + v);
					x1.imag(0.0);

					x2.real(-b / (3.0 * a) - (u + v) / 2.0);
					x2.imag(std::sqrt(3.0)*(u - v) / 2.0);

					x3.real(-b / (3.0 * a) - (u + v) / 2.0);
					x3.imag(-std::sqrt(3.0)*(u - v) / 2.0);
					return { x1, x2, x3 };
				}
				else if (delta == 0) {
					auto u = real_cubic_root(-q / 2.0 + std::pow(delta, 0.5));
					auto v = real_cubic_root(-q / 2.0 - std::pow(delta, 0.5));
					std::complex<double> x1, x2, x3;
					x1.real(-b / (3.0 * a) + u + v);
					x1.imag(0.0);

					x2.real(-b / (3.0 * a) - (u + v) / 2.0);
					x2.imag(0.0);

					x3.real(-b / (3.0 * a) - (u + v) / 2.0);
					x3.imag(0.0);
					return { x1, x2, x3 };
				}
				else {
					auto phi = std::acos(-q / (2.0 * std::sqrt(std::pow(-p / 3.0, 3.0))));
					std::complex<double> x1, x2, x3;
					x1.real(-b / (3.0 * a) + 2.0 * std::sqrt(-p / 3.0)*std::cos(phi / 3.0));
					x1.imag(0.0);

					x2.real(-b / (3.0 * a) + 2.0 * std::sqrt(-p / 3.0)*std::cos((phi + 2.0 * M_PI) / 3.0));
					x2.imag(0.0);

					x3.real(-b / (3.0 * a) + 2 * std::sqrt(-p / 3.0)*std::cos((phi - 2.0 * M_PI) / 3.0));
					x3.imag(0.0);
					return { x1, x2, x3 };
				}
			};

			auto quad_roots = [cubic_roots](double a, double b, double c, double d, double e)->std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>>
			{
				auto delta = 256.0 * std::pow(a*e, 3.0) - 192.0 * std::pow(a*e, 2.0) * b * d - 128.0 * std::pow(a*c*e, 2.0) + 144.0 * std::pow(a*d, 2.0) * c*e - 27.0 * std::pow(a*std::pow(d, 2.0), 2.0) + 144.0 * std::pow(b*e, 2.0) * a*c - 6.0 * std::pow(b*d, 2.0) * a*e - 80.0 * a*b*std::pow(c, 2.0) * d*e + 18.0 * a*b*c*std::pow(d, 3.0) + 16.0 * a*std::pow(c, 4.0) * e - 4.0 * a*std::pow(c, 3.0) * std::pow(d, 2.0) - 27.0 * std::pow(b, 4.0) * std::pow(e, 2.0) + 18.0 * std::pow(b, 3.0) * c*d*e - 4.0 * std::pow(b*d, 3.0) - 4.0 * std::pow(b, 2.0) * std::pow(c, 3.0) * e + std::pow(b*c*d, 2.0);
				auto P = 8.0 * a*c - 3.0 * std::pow(b, 2.0);
				auto R = std::pow(b, 3.0) + 8.0 * d*std::pow(a, 2.0) - 4 * a*b*c;
				auto Delta0 = std::pow(c, 2.0) - 3.0 * b*d + 12.0 * a*e;
				auto D = 64.0 * std::pow(a, 3.0) * e - 16.0 * std::pow(a*c, 2.0) + 16.0 * a*std::pow(b, 2.0) * c - 16 * std::pow(a, 2.0) * b*d - 3 * std::pow(b, 4.0);

				auto A = c / a - (3.0 / 8.0)*std::pow(b / a, 2.0);
				auto B = d / a - b * c / (2.0 * std::pow(a, 2.0)) + 1.0 / 8.0 * std::pow(b / a, 3.0);
				auto C = e / a - b * d / (4.0 * std::pow(a, 2.0)) + std::pow(b, 2.0) * c / (16.0 * std::pow(a, 3.0)) - 3.0 / 256.0 * std::pow(b / a, 4.0);

				auto coef1 = 8.0;
				auto coef2 = -4.0 * A;
				auto coef3 = -8.0 * C;
				auto coef4 = 4.0 * A*C - std::pow(B, 2.0);


				std::tuple<std::complex<double>, std::complex<double>, std::complex<double>> cub_roots = cubic_roots(coef1, coef2, coef3, coef4);

				auto x10 = 2.0 * std::get<0>(cub_roots) - A;
				auto x11 = -2.0 * std::get<0>(cub_roots) - A + 2.0 * B / (std::sqrt(2.0 * std::get<0>(cub_roots) - A));
				auto x12 = std::sqrt(2.0 * std::get<0>(cub_roots) - A);
				std::complex<double> x = std::sqrt(std::complex<double>(-1.0));

				std::complex<double> y1 = -(1.0 / 2.0)*std::pow(2.0 * std::get<0>(cub_roots) - A, 0.5) + (1.0 / 2.0)*std::pow(-2.0 * std::get<0>(cub_roots) - A + 2.0 * B / (std::sqrt(2.0 * std::get<0>(cub_roots) - A)), 0.5);
				std::complex<double> y2 = -(1.0 / 2.0)*std::pow(2.0 * std::get<0>(cub_roots) - A, 0.5) - (1.0 / 2.0)*std::pow(-2.0 * std::get<0>(cub_roots) - A + 2.0 * B / (std::sqrt(2.0 * std::get<0>(cub_roots) - A)), 0.5);
				std::complex<double> y3 = (1.0 / 2.0)*std::pow(2.0 * std::get<0>(cub_roots) - A, 0.5) + (1.0 / 2.0)*std::pow(-2.0 * std::get<0>(cub_roots) - A - 2.0 * B / (std::sqrt(2.0 * std::get<0>(cub_roots) - A)), 0.5);
				std::complex<double> y4 = (1.0 / 2.0)*std::pow(2.0 * std::get<0>(cub_roots) - A, 0.5) - (1.0 / 2.0)*std::pow(-2.0 * std::get<0>(cub_roots) - A - 2.0 * B / (std::sqrt(2.0 * std::get<0>(cub_roots) - A)), 0.5);

				std::complex<double> x1 = y1 - b / (4.0 * a);
				std::complex<double> x2 = y2 - b / (4.0 * a);
				std::complex<double> x3 = y3 - b / (4.0 * a);
				std::complex<double> x4 = y4 - b / (4.0 * a);

				return { x1, x2, x3, x4 };

			};

			auto clean_roots = [](std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>> roots, double b) -> std::vector<double>
			{
				std::vector<double> cleanroots;
				if (std::get<0>(roots).imag() == 0) {
					if (std::get<0>(roots).real() > b) {
						cleanroots.push_back(std::get<0>(roots).real());
					}
				}
				if (std::get<1>(roots).imag() == 0) {
					if (std::get<1>(roots).real() > b) {
						cleanroots.push_back(std::get<1>(roots).real());
					}
				}
				if (std::get<2>(roots).imag() == 0) {
					if (std::get<2>(roots).real() > b) {
						cleanroots.push_back(std::get<2>(roots).real());
					}
				}
				if (std::get<3>(roots).imag() == 0) {
					if (std::get<3>(roots).real() > b) {
						cleanroots.push_back(std::get<3>(roots).real());
					}
				}
				return cleanroots;
			};

			auto k1 = epsilon * bmix.val;
			auto k2 = sigma * bmix.val;
			auto AP = R * T.val;
			auto BP = 2.0 * R*T.val*(k2 + k1) - 2.0 * amix.val;
			auto CP = R * T.val * (std::pow(k2, 2.0) + 4.0 * k1 * k2 + std::pow(k1, 2.0)) - amix.val * (k1 + k2 - 4.0 * bmix.val);
			auto DP = 2 * R*T.val*(k1*std::pow(k2, 2.0) + std::pow(k1, 2.0) * k2) - 2.0 * amix.val * (std::pow(bmix.val, 2.0) - k2 * bmix.val - k1 * bmix.val);
			auto EP = R * T.val*std::pow(k1, 2.0) * std::pow(k2, 2.0) - amix.val * (k1 + k2)*std::pow(bmix.val, 2.0);

			auto ro = quad_roots(AP, BP, CP, DP, EP);

			auto user_roots = clean_roots(ro, bmix.val);

			if (user_roots.size() != 2) {
				should_be_gas_gibbs = true;
			}
			else {
				auto pmin = p(user_roots[0]);
				auto pmax = p(user_roots[1]);

				if (pmin > pmax) {
					auto temp = pmax;
					pmax = pmin;
					pmin = temp;
				}

				if (P.val > pmin && P.val < pmax) {
					std::cout << "indeterminado" << std::endl;
					should_be_gas_gibbs = true;
				}
				if (P.val < pmin) {
					should_be_gas_gibbs = true;
				}
				if (P.val > pmax) {
					should_be_gas_gibbs = false;
				}
			}
		}

		if ((real_roots.size() != 1)){// && minimaze_gibbs_energy) {
			std::array<ChemicalScalar, 2> Gs;
			std::vector<ChemicalVector> ln_phis;
			ln_phis.push_back(ChemicalVector(2));
			ln_phis.push_back(ChemicalVector(2));
			std::array<ChemicalScalar, 2> Gibbs_energy;
			for (auto i = 0; i < 2; i++) {
				// Calculate the partial derivatives of each Z (dZdT, dZdP, dZdn)
				const double factor = -1.0 / (3 * Zs[i].val*Zs[i].val + 2 * A.val*Zs[i].val + B.val);
				Zs[i].ddT = factor * (A.ddT*Zs[i].val*Zs[i].val + B.ddT*Zs[i].val + C.ddT);
				Zs[i].ddP = factor * (A.ddP*Zs[i].val*Zs[i].val + B.ddP*Zs[i].val + C.ddP);
				for (unsigned j = 0; j < nspecies; ++j)
					Zs[i].ddn[j] = factor * (A.ddn[j] * Zs[i].val*Zs[i].val + B.ddn[j] * Zs[i].val + C.ddn[j]);

				// Calculate the partial temperature derivative of Z
				const ChemicalScalar ZT = -(AT*Zs[i] * Zs[i] + BT * Zs[i] + CT) / (3 * Zs[i] * Zs[i] + 2 * A*Zs[i] + B);

				// Calculate the integration factor I and its temperature derivative IT
				ChemicalScalar I;
				if (epsilon != sigma) I = log((Zs[i] + sigma * beta) / (Zs[i] + epsilon * beta)) / (sigma - epsilon);
				else I = beta / (Zs[i] + epsilon * beta);

				// Calculate the temperature derivative IT of the integration factor I
				ChemicalScalar IT;
				if (epsilon != sigma) IT = ((ZT + sigma * betaT) / (Zs[i] + sigma * beta) - (ZT + epsilon * betaT) / (Zs[i] + epsilon * beta)) / (sigma - epsilon);
				else IT = I * (betaT / beta - (ZT + epsilon * betaT) / (Zs[i] + epsilon * beta));

				Gs[i] = (R*T*(Zs[i] - 1 - log(Zs[i] - beta) - q * I));
				//ln_phis[i] = Gs[i] / (R*T);
			}
			/*
			for (auto i = 0; i < 2; i++) {
				for (auto j = 0; j < x.val.size(); j++) {
					Gibbs_energy[i].val += ln_phis[i][j].val * x[j].val;
				}
			}
			if (Gibbs_energy[0] > Gibbs_energy[1]) {
				Z.val = *std::max_element(std::begin(real_roots), std::end(real_roots));
			}
			else {
				Z.val = *std::min_element(std::begin(real_roots), std::end(real_roots));
			}
			*/
			if (Gs[0] > Gs[1]) {
				should_be_gas_gibbs = true;
			}
			else {
				should_be_gas_gibbs = false;
			}
		}


		/*
		if (isvapor && !should_be_gas ) {
			Z = Z * 1000;
		}

		if (!isvapor && should_be_gas) {
			Z = Z / 1000;
		}
		*/

        /*
        double greatest_root = (modulo1 > modulo2) ? modulo1 : modulo2;
        double smallest_root = (modulo1 > modulo2) ? modulo2 : modulo1;
        */
        /*
        double greatest_root = (first_root.real() > second_root.real()) ? first_root.real() : second_root.real();
        double smallest_root = (first_root.real() > second_root.real()) ? second_root.real() : first_root.real();
        */

        //std::cout << "isvapor: " << isvapor << "\ngreatest_root: " << *std::max_element(std::begin(real_roots), std::end(real_roots)) << "\nmodulo_smallest_root: " << *std::min_element(std::begin(real_roots), std::end(real_roots)) << std::endl;
        
        //std::cout << "used" << Z.val <<std::endl;
        
        

        //std::cout << "old newton: " << newton(f, Z0, tolerance, maxiter) << std::endl;

        //Z.val = isvapor ? greatest_root : smallest_root;


        //Z.val = isvapor ? bounds.first.real() : bounds.second.real();

       

        /*
        if ((greatest_root == smallest_root) && (Z.val > 0.9)) {
            if (!isvapor) {
                Z.val = greatest_root + 10;
            }
        }
        if ((greatest_root == smallest_root) && (Z.val < 0.3)) {
            if (isvapor) {
                Z.val = greatest_root - 10;
            }
        }
        */
        // Calculate the partial derivatives of Z (dZdT, dZdP, dZdn)
        const double factor = -1.0/(3*Z.val*Z.val + 2*A.val*Z.val + B.val);
        Z.ddT = factor * (A.ddT*Z.val*Z.val + B.ddT*Z.val + C.ddT);
        Z.ddP = factor * (A.ddP*Z.val*Z.val + B.ddP*Z.val + C.ddP);
        for(unsigned i = 0; i < nspecies; ++i)
            Z.ddn[i] = factor * (A.ddn[i]*Z.val*Z.val + B.ddn[i]*Z.val + C.ddn[i]);

        // Calculate the partial temperature derivative of Z
        const ChemicalScalar ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B);

        // Calculate the integration factor I and its temperature derivative IT
        ChemicalScalar I;
        if(epsilon != sigma) I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon);
                        else I = beta/(Z + epsilon*beta);

        // Calculate the temperature derivative IT of the integration factor I
        ChemicalScalar IT;
        if(epsilon != sigma) IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon);
                        else IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta));


        
        //std::cout << "isvapor: " << isvapor << " fugacity" << exp(ln_phi) << std::endl;

    


		/*
		const auto dpdv = [&](double V) -> std::tuple<double, double>
		{
			const double val = ((amix.val*(2.0*bmix.val+2.0*V))/(std::pow(-(bmix.val*bmix.val)+2.0*bmix.val*V+V*V,2.0))) - ((R*T.val)/std::pow(V-bmix.val,2.0));
			const double g1 = (2.0 * R*T.val) / (std::pow(V - bmix.val, 3.0));
			const double g2 = (2.0 * (std::pow(2.0 * bmix.val + 2 * V, 2) ) )/ std::pow(-(bmix.val*bmix.val) + 2 * bmix.val*V + V * V, 3);
			const double g3 = 2.0 / std::pow(-(bmix.val*bmix.val) + 2.0 * bmix.val*V + V * V, 2.0);
			const double grad = g1 - amix.val*(g2 - g3);
			return std::make_tuple(val, grad);
		};
		*/
		/*
		auto Vmin = 1.0e-5;// newton(dpdv, 1.0e-5, tolerance, maxiter);
		auto Vmax = 0.00025;// newton(dpdv, 0.00025, tolerance, maxiter);

		
		std::cout << pmin << " " << pmax << std::endl;
		std::cout << Vmin << " " << Vmax << std::endl;
		std::cout << P.val << std::endl;
		std::cout << V.val << std::endl;

		auto teste1 = Gs[0].val;
		auto teste2 = Gs[1].val;
		

		*/
		
		//bool should_be_gas = workanalisys::should_be_gas(V_gas, V_liq, amix.val, bmix.val, R, T.val, P.val);

		

		


        // Calculate the partial molar Zi for each species
		V = Z * R*T / P;
        G_res = R*T*(Z - 1 - log(Z - beta) - q*I);
        H_res = R*T*(Z - 1 + T*qT*I);
        Cp_res = R*T*(ZT + qT*I + T*qTT + T*qT*IT) + H_res/T;

        const ChemicalScalar dPdT = P*(1.0/T + ZT/Z);
        const ChemicalScalar dVdT = V*(1.0/T + ZT/Z);

        Cv_res = Cp_res - T*dPdT*dVdT + R;
		

        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double bi = bbar[i];
            const ThermoScalar betai = P*bi/(R*T);
            const ChemicalScalar ai = abar[i];
            const ChemicalScalar aiT = abarT[i];
            const ChemicalScalar qi = q*(1 + ai/amix - bi/bmix);
            const ChemicalScalar qiT = qi*qT/q + q*(aiT - ai*amixT/amix)/amix;
            const ThermoScalar Ai = (epsilon + sigma - 1.0)*betai - 1.0;
            const ChemicalScalar Bi = (epsilon*sigma - epsilon - sigma)*(2*beta*betai - beta*beta) - (epsilon + sigma - q)*(betai - beta) - (epsilon + sigma - qi)*beta;
            const ChemicalScalar Ci = -3*sigma*epsilon*beta*beta*betai + 2*epsilon*sigma*beta*beta*beta - (epsilon*sigma + qi)*beta*beta - 2*(epsilon*sigma + q)*(beta*betai - beta*beta);
            const ChemicalScalar Zi = -(Ai*Z*Z + (Bi + B)*Z + Ci + 2*C)/(3*Z*Z + 2*A*Z + B);

			ChemicalScalar Ii;
            if(epsilon != sigma) Ii = I + ((Zi + sigma*betai)/(Z + sigma*beta) - (Zi + epsilon*betai)/(Z + epsilon*beta))/(sigma - epsilon);
                            else Ii = I * (1 + betai/beta - (Zi + epsilon*betai)/(Z + epsilon*beta));

            Vi[i] = R*T*Zi/P;
            Gi_res[i] = R*T*(Zi - (Zi - betai)/(Z - beta) - log(Z - beta) - qi*I - q*Ii + q*I);
            Hi_res[i] = R*T*(Zi - 1 + T*(qiT*I + qT*Ii - qT*I));
            ln_phi[i] = Gi_res[i]/(R*T);
        }

        //test if phase is stable - applaying Volume Methods (Bennett, Jim, and Kurt AG Schmidt. "Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations." Energy & Fuels 31.4 (2016): 3370-3379.)
		
		//bool should_be_gas = VolumeMethod::should_be_gas(V.val, bmix.val);

		//bool should_be_gas = CriticalPointMethods::should_be_gas(amix.val, bmix.val, R, Omega, Psi, Z.val, V.val, T.val);
		
		//bool should_be_gas = IsothermalCompressibilityMethods::should_be_gas(T, P, Z, R);
		
		

		
		switch (phaseidMethod)
		{
		case phaseIdentificationMethod::workanalisys:
			should_be_gas = workanalisys::should_be_gas(V_gas.val, V_liq.val, amix.val, bmix.val, R, T.val, P.val);
			break;
		case phaseIdentificationMethod::IsothermalCompressibilityMethods:
			should_be_gas = IsothermalCompressibilityMethods::should_be_gas(T, P, Z, R);
			break;
		case phaseIdentificationMethod::CriticalPointMethods:
			should_be_gas = CriticalPointMethods::should_be_gas(amix.val, bmix.val, R, Omega, Psi, Z.val, V.val, T.val);
			break;
		case phaseIdentificationMethod::VolumeMethod:
			should_be_gas = VolumeMethod::should_be_gas(V.val, bmix.val);
			break;
		case phaseIdentificationMethod::workanalisysPengRobinson:
			should_be_gas = workanalisys_PengRobinson::should_be_gas(V_gas.val, V_liq.val, amix.val, bmix.val, R, T.val, P.val);
			break;
		case phaseIdentificationMethod::Gibbs_residual_based:
			should_be_gas = should_be_gas_gibbs;
			break;
		default:
			should_be_gas = isvapor;
			break;
		}

		auto should_be_gas1 = workanalisys::should_be_gas(V_gas.val, V_liq.val, amix.val, bmix.val, R, T.val, P.val);
		auto should_be_gas2 = IsothermalCompressibilityMethods::should_be_gas(T, P, Z, R);
		auto should_be_gas3 = CriticalPointMethods::should_be_gas(amix.val, bmix.val, R, Omega, Psi, Z.val, V.val, T.val);
		auto should_be_gas4 = VolumeMethod::should_be_gas(V.val, bmix.val);
		auto should_be_gas5 = workanalisys_PengRobinson::should_be_gas(V_gas.val, V_liq.val, amix.val, bmix.val, R, T.val, P.val);

		/*
		if (should_be_gas4 != should_be_gas_gibbs) {
			std::cout << "deferentes definições de estado" << std::endl;
			std::cout << (isvapor ? "gas" : "oil") << std::endl;
			std::cout << "volume method clasified as: " << (should_be_gas4 ? "gas" : "oil") << std::endl;
			std::cout << "gibbs method clasified as: " << (should_be_gas_gibbs ? "gas" : "oil") << std::endl;
			std::cout << "quantity = " << quantity << std::endl;
		}
		*/
		if (!isvapor) {
			quantity--;
		}
		//std::cout << "QUANTITYYYY" << quantity << std::endl;
		
        if (isvapor && !should_be_gas && (quantity > 0) ){// && (quantity > 390)) {
            //std::cout << "zerando gas" << std::endl;
			result.molar_volume = 0.0;
            result.residual_molar_gibbs_energy = 0.0;
            result.residual_molar_enthalpy = 0.0;
            result.residual_molar_heat_capacity_cp = 0.0;
            result.residual_molar_heat_capacity_cv = 0.0;
            result.partial_molar_volumes.fill(0.0);
            result.residual_partial_molar_gibbs_energies.fill(0.0);
            result.residual_partial_molar_enthalpies.fill(0.0);
			//result.ln_fugacity_coefficients.fill(0.0, 1);// , 0.020634908961432240, -1.4948782256645683e-06, +9.e30);//result.ln_fugacity_coefficients.fill(100);						
			result.ln_fugacity_coefficients.fill(100.0);
        }
        
        if (!isvapor && should_be_gas && (quantity > 0)){// && (quantity > 190)) {
            //std::cout << "zerando oil" << std::endl;
			result.molar_volume = 0.0;
            result.residual_molar_gibbs_energy = 0.0;
            result.residual_molar_enthalpy = 0.0;
            result.residual_molar_heat_capacity_cp = 0.0;
            result.residual_molar_heat_capacity_cv = 0.0;
            result.partial_molar_volumes.fill(0.0);
            result.residual_partial_molar_gibbs_energies.fill(0.0);
            result.residual_partial_molar_enthalpies.fill(0.0);
			//result.ln_fugacity_coefficients.fill(0.0, 1);//, 0.020634908961432240, -1.4948782256645683e-06, +9.e30);
			result.ln_fugacity_coefficients.fill(100.0);
        }
		

//		std::cout << "is " << (isvapor ? "vapor" : "oil") << std::endl;
//		std::cout << (should_be_gas ? "should be gas" : "shouldn't be gas") << std::endl;

        return result;
    }
};

CubicEOS::Result::Result()
{}

CubicEOS::Result::Result(unsigned nspecies)
: molar_volume(nspecies),
  residual_molar_gibbs_energy(nspecies),
  residual_molar_enthalpy(nspecies),
  residual_molar_heat_capacity_cp(nspecies),
  residual_molar_heat_capacity_cv(nspecies),
  partial_molar_volumes(nspecies),
  residual_partial_molar_gibbs_energies(nspecies),
  residual_partial_molar_enthalpies(nspecies),
  ln_fugacity_coefficients(nspecies)
{}

CubicEOS::CubicEOS(unsigned nspecies)
: pimpl(new Impl(nspecies))
{}

CubicEOS::CubicEOS(const CubicEOS& other)
: pimpl(new Impl(*other.pimpl))
{}

CubicEOS::~CubicEOS()
{}

auto CubicEOS::operator=(CubicEOS other) -> CubicEOS&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto CubicEOS::numSpecies() const -> unsigned
{
    return pimpl->nspecies;
}

auto CubicEOS::setModel(Model model) -> void
{
    pimpl->model = model;
}

/*
auto CubicEOS::setPhaseId(phaseidentification id) -> void
{
	pimpl->phaseid = id;
}
*/

auto CubicEOS::setPhaseAsLiquid() -> void
{
    pimpl->isvapor = false;
}

auto CubicEOS::setPhaseAsVapor() -> void
{
    pimpl->isvapor = true;
}

auto CubicEOS::setMinimazeGibbsEnergytrue() -> void
{
	pimpl->minimaze_gibbs_energy = true;
}

auto CubicEOS::setMinimazeGibbsEnergyfalse() -> void
{
	pimpl->minimaze_gibbs_energy = false;
}

auto CubicEOS::setCriticalTemperatures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "temperatures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical temperatures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "temperatures of the gases.");

    pimpl->critical_temperatures = values;
}

auto CubicEOS::setCriticalPressures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "pressures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical pressures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "pressures of the gases.");

    pimpl->critical_pressures = values;
}

auto CubicEOS::setAcentricFactors(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the acentric "
        "factors of the species in CubicEOS.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " values were given.");

    pimpl->acentric_factors = values;
}

auto CubicEOS::setInteractionParamsFunction(const InteractionParamsFunction& func) -> void
{
    pimpl->calculate_interaction_params = func;
}

auto CubicEOS::operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> Result
{
    return pimpl->operator()(T, P, x);
}

} // namespace Reaktoro
