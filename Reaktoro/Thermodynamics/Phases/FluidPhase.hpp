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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/EOS/CubicEOS.hpp>

namespace Reaktoro {

// Forward declarations
class FluidMixture;

/// Class that defines a fluid (gaseous or liquid) phase
class FluidPhase : public Phase
{
public:
    ///Construct an FluidPhase instance with given name and PhaseType
    FluidPhase(const std::string& name, PhaseType type);

    /// Construct an FluidPhase instance with given fluid (gaseous or liquid) mixture, name and type.
    /// The Peng-Robinson equation of state is chosen by default to calculate the
    /// thermodynamic and chemical properties of this FluidPhase object.
    explicit FluidPhase(const FluidMixture& mixture, const std::string& name, PhaseType type);

    /// Set the chemical model of the phase with the ideal gas equation of state.
    /// This model only supports a gaseous phase. Using it in a FluidPhase that is not a
    /// PhaseType::Gas will result in a runtime error.
    auto setChemicalModelIdeal() -> FluidPhase&;

    /// Set the chemical model of the phase with the van der Waals equation of state.
    ///
    /// Reference: *van der Waals, J.D. (1910). Nobel Lectures in Physics. pp. 254-265*.
    auto setChemicalModelVanDerWaals() -> FluidPhase&;

    /// Set the chemical model of the phase with the Redlich-Kwong equation of state.
    ///
    /// Reference: *Redlich, O., Kwong, J.N.S. (1949). On The Thermodynamics of Solutions. Chem. Rev. 44(1) 233–244*.
    auto setChemicalModelRedlichKwong() -> FluidPhase&;

    /// Set the chemical model of the phase with the Soave-Redlich-Kwong equation of state.
    ///
    /// Reference: *Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci., 27, 1197-1203*.
    auto setChemicalModelSoaveRedlichKwong() -> FluidPhase&;

    /// Set the chemical model of the phase with the Peng-Robinson equation of state.
    ///
    /// Reference: *Peng, D.Y., Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial and Engineering Chemistry: Fundamentals 15: 59–64*.
    auto setChemicalModelPengRobinson(CubicEOS::Params params = {}) -> FluidPhase&;

    /// Set the chemical model of the phase with a Cubic equation of state. The specific type
    /// can be passed as a parameter (the default in Peng-Robinson).
    ///
    /// References:
    ///     *Peng, D.Y., Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial and Engineering Chemistry: Fundamentals 15: 59-64*.
    ///     *Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci., 27, 1197-1203*.
    ///     *Redlich, O., Kwong, J.N.S. (1949). On The Thermodynamics of Solutions. Chem. Rev. 44(1) 233-244*.
    ///     *van der Waals, J.D. (1910). Nobel Lectures in Physics. pp. 254-265*.
    auto setChemicalModelCubicEOS(CubicEOS::Params params = {}) -> FluidPhase&;

    /// Set the chemical model of the phase with the Spycher et al. (2003) equation of state.
    /// This model only supports the gaseous species `H2O(g)` and `CO2(g)`. Any other species
    /// will result in a runtime error.
    ///
    /// Reference: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
    /// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
    /// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
    auto setChemicalModelSpycherPruessEnnis() -> FluidPhase&;

    /// Set the chemical model of the phase with the Spycher and Reed (1988) equation of state.
    /// This model only supports the gaseous species `H2O(g)`, `CO2(g)`, and CH4(g). Any other
    /// species will result in a runtime error.
    ///
    /// Reference: *Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2,
    /// CH4, H2O and of H2O--CO2--CH4 mixtures: A virial equation treatment for
    /// moderate pressures and temperatures applicable to calculations of
    /// hydrothermal boiling. Geochimica et Cosmochimica Acta, 52(3), 739–749*.
    auto setChemicalModelSpycherReed() -> FluidPhase&;

    /// Return a const reference of the FluidMixture instance
    auto mixture() const -> const FluidMixture&;

    /// Return the FluidMixture instance
    auto mixture()->FluidMixture&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
