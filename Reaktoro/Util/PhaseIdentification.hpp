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

#pragma once

#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ChemicalScalar.hpp>

namespace Reaktoro {

    // Forward declarations
    enum class PhaseType;

namespace PhaseIdentification {
    
    /// Return a PhaseType that says if the phase is a Liquid or Gas based on Volume Method
    /// @param Volume Phase volume
    /// @param b coefficient of EoS
    /// @return The type of the phase
    /// 
    /// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
    auto volumeMethod(ChemicalScalar& Volume, ChemicalScalar& b)->PhaseType;

    /// Return a PhaseType that says if the phase is a Liquid or Gas based on Isothermal Compressibility
    /// @param Temperature Phase temperature
    /// @param Pressure Phase pressure
    /// @param Z Phase compressibility 
    /// @return The type of the phase
    /// 
    /// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
    auto isothermalCompressibilityMethod(const ThermoScalar& Temperature, const ThermoScalar& Pressure, const ChemicalScalar& Z)->PhaseType;
    
    /// Return a PhaseType that says if the phase is a Liquid or Gas based on gibbs residual energy and equation of state
    /// if size(Z) > 2, the values od Gibbs residual energy is compared. It is a liquid phase if Gibbs residual energy of Z_min is the smallest and gaseous if 
    /// Gibbs residual energy of Z_max is the smallest.
    /// if size(Z) == 1, the pressue is compared with the local P_min e local P_max of the EoS. It is a liquid phase if P \textgreater P_min and gaseous if P \textless Pmax
    /// @param Pressure Phase pressure
    /// @param Temperature Phase temperature
    /// @param amix attractive parameter
    /// @param Z Phase compressibility 
    /// @return The type of the phase
    /// 
    /// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
    auto gibbsEnergyAndEquationOfStateMethod(const ThermoScalar& Pressure, const ThermoScalar& Temperature, const ChemicalScalar& amix, 
                                             const ChemicalScalar& bmix, const ChemicalScalar& A, const ChemicalScalar& B, const ChemicalScalar& C, 
                                             std::vector<ChemicalScalar> Z, const double epsilon, const double sigma)->PhaseType;


}

}
