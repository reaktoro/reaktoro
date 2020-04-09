// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The primary chemical property data of the phases and species in a chemical system.
/// @see ChemicalProps
struct ChemicalPropsData
{
    /// The temperature of the system (in K).
    real T = {};

    /// The pressure of the system (in Pa).
    real P = {};

    /// The amounts of each species in the system (in mol).
    ArrayXr n;

    /// The sum of species amounts in each phase of the system (in mol).
    ArrayXr amount;

    /// The mole fractions of the species in the system (in mol/mol).
    ArrayXr x;

    /// The standard molar Gibbs energies of the species in the system (in J/mol).
    ArrayXr G0;

    /// The standard molar enthalpies of the species in the system (in J/mol).
    ArrayXr H0;

    /// The standard molar volumes of the species in the system (in m3/mol).
    ArrayXr V0;

    /// The standard molar isobaric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cp0;

    /// The standard molar isochoric heat capacities of the species in the system (in J/(mol·K)).
    ArrayXr Cv0;

    /// The excess molar volume of each phase in the system (in m3/mol).
    ArrayXr Vex;

    /// The temperature derivative at constant pressure of the excess molar volume of each phase in the system (in m3/(mol*K)).
    ArrayXr VexT;

    /// The pressure derivative at constant temperature of the excess molar volume of each phase in the system (in m3/(mol*Pa)).
    ArrayXr VexP;

    /// The excess molar Gibbs energy of each phase in the system (in units of J/mol).
    ArrayXr Gex;

    /// The excess molar enthalpy of each phase in the system (in units of J/mol).
    ArrayXr Hex;

    /// The excess molar isobaric heat capacity of each phase in the system (in units of J/(mol*K)).
    ArrayXr Cpex;

    /// The excess molar isochoric heat capacity of each phase in the system (in units of J/(mol*K)).
    ArrayXr Cvex;

    /// The activity coefficients (natural log) of the species in the system.
    ArrayXr ln_g;

    /// The activities (natural log) of the species in the system.
    ArrayXr ln_a;
};

/// The chemical properties of the phases and species in a chemical system.
class ChemicalProps
{
public:
    /// Construct a ChemicalProps object.
    explicit ChemicalProps(const ChemicalSystem& system);

    /// Construct a ChemicalProps object.
    ChemicalProps(const ChemicalSystem& system, const ChemicalPropsData& data);

    /// Return the chemical system associated with these chemical properties.
    auto system() const -> const ChemicalSystem&;

    /// Return the primary chemical property data of the system from which others are calculated.
    auto data() const -> const ChemicalPropsData&;

    /// Return the chemical properties of a phase with given index.
    auto phase(Index idx) const -> ChemicalPropsPhaseConstRef;

    /// Return the chemical properties of a phase with given index.
    auto phase(Index idx) -> ChemicalPropsPhaseRef;

    /// Return the temperature of the system (in K).
    auto temperature() const -> real;

    /// Return the pressure of the system (in Pa).
    auto pressure() const -> real;

    /// Return the amounts of each species in the system (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the mole fractions of the species in the system.
    auto moleFractions() const -> ArrayXrConstRef;

    /// Return the ln activity coefficients of the species in the system.
    auto lnActivityCoefficients() const -> ArrayXrConstRef;

    /// Return the ln activities of the species in the system.
    auto lnActivities() const -> ArrayXrConstRef;

    /// Return the chemical potentials of the species in the system (in J/mol).
    auto chemicalPotentials() const -> ArrayXr;

    /// Return the standard partial molar Gibbs energies of the species in the system (in J/mol).
    auto standardGibbsEnergies() const -> ArrayXrConstRef;

    /// Return the standard partial molar enthalpies of the species in the system (in J/mol).
    auto standardEnthalpies() const -> ArrayXrConstRef;

    /// Return the standard partial molar volumes of the species in the system (in m3/mol).
    auto standardVolumes() const -> ArrayXrConstRef;

    /// Return the standard partial molar entropies of the species in the system (in J/(mol*K)).
    auto standardEntropies() const -> ArrayXr;

    /// Return the standard partial molar internal energies of the species in the system (in J/mol).
    auto standardInternalEnergies() const -> ArrayXr;

    /// Return the standard partial molar Helmholtz energies of the species in the system (in J/mol).
    auto standardHelmholtzEnergies() const -> ArrayXr;

    /// Return the standard partial molar isobaric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const -> ArrayXrConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const -> ArrayXrConstRef;

private:
    /// The chemical system associated with these chemical properties.
    ChemicalSystem sys;

    /// The primary chemical property data of the system from which others are calculated.
    ChemicalPropsData props;
};

} // namespace Reaktoro
