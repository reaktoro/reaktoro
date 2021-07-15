// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "Prop.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

namespace Reaktoro {

Prop::Prop(const PropFn& propfn)
: propfn(propfn)
{}

auto Prop::eval(const ChemicalProps& props) -> real
{
    return propfn(props);
}

auto Prop::operator()(const ChemicalProps& props) -> real
{
    return propfn(props);
}

auto Prop::elementAmount(const ChemicalSystem& system, const String& element) -> Prop
{
    errorif(true, "Prop::elementAmount has not been implemented yet."); // TODO: Implement method Prop::elementAmount.
}

auto Prop::elementAmountInPhase(const ChemicalSystem& system, const String& element, const String& phase) -> Prop
{
    errorif(true, "Prop::elementAmountInPhase has not been implemented yet."); // TODO: Implement method Prop::elementAmountInPhase.
}

auto Prop::elementMass(const ChemicalSystem& system, const String& element) -> Prop
{
    errorif(true, "Prop::elementMass has not been implemented yet."); // TODO: Implement method Prop::elementMass.
}

auto Prop::elementMassInPhase(const ChemicalSystem& system, const String& element, const String& phase) -> Prop
{
    errorif(true, "Prop::elementMassInPhase has not been implemented yet."); // TODO: Implement method Prop::elementMassInPhase.
}

auto Prop::speciesAmount(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::speciesAmount has not been implemented yet."); // TODO: Implement method Prop::speciesAmount.
}

auto Prop::speciesMass(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::speciesMass has not been implemented yet."); // TODO: Implement method Prop::speciesMass.
}

auto Prop::moleFraction(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::moleFraction has not been implemented yet."); // TODO: Implement method Prop::moleFraction.
}

auto Prop::activityCoefficient(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::activityCoefficient has not been implemented yet."); // TODO: Implement method Prop::activityCoefficient.
}

auto Prop::activity(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::activity has not been implemented yet."); // TODO: Implement method Prop::activity.
}

auto Prop::chemicalPotential(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::chemicalPotential has not been implemented yet."); // TODO: Implement method Prop::chemicalPotential.
}

auto Prop::standardVolume(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardVolume has not been implemented yet."); // TODO: Implement method Prop::standardVolume.
}

auto Prop::standardGibbsEnergy(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardGibbsEnergy has not been implemented yet."); // TODO: Implement method Prop::standardGibbsEnergy.
}

auto Prop::standardEnthalpy(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardEnthalpy has not been implemented yet."); // TODO: Implement method Prop::standardEnthalpy.
}

auto Prop::standardEntropy(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardEntropy has not been implemented yet."); // TODO: Implement method Prop::standardEntropy.
}

auto Prop::standardInternalEnergy(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardInternalEnergy has not been implemented yet."); // TODO: Implement method Prop::standardInternalEnergy.
}

auto Prop::standardHelmholtzEnergy(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardHelmholtzEnergy has not been implemented yet."); // TODO: Implement method Prop::standardHelmholtzEnergy.
}

auto Prop::standardHeatCapacitiesConstP(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardHeatCapacitiesConstP has not been implemented yet."); // TODO: Implement method Prop::standardHeatCapacitiesConstP.
}

auto Prop::standardHeatCapacitiesConstV(const ChemicalSystem& system, const String& species) -> Prop
{
    errorif(true, "Prop::standardHeatCapacitiesConstV has not been implemented yet."); // TODO: Implement method Prop::standardHeatCapacitiesConstV.
}

auto Prop::phaseAmount(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseAmount has not been implemented yet."); // TODO: Implement method Prop::phaseAmount.
}

auto Prop::phaseMass(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseMass has not been implemented yet."); // TODO: Implement method Prop::phaseMass.
}

auto Prop::phaseVolume(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseVolume has not been implemented yet."); // TODO: Implement method Prop::phaseVolume.
}

auto Prop::phaseGibbsEnergy(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseGibbsEnergy has not been implemented yet."); // TODO: Implement method Prop::phaseGibbsEnergy.
}

auto Prop::phaseEnthalpy(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseEnthalpy has not been implemented yet."); // TODO: Implement method Prop::phaseEnthalpy.
}

auto Prop::phaseEntropy(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseEntropy has not been implemented yet."); // TODO: Implement method Prop::phaseEntropy.
}

auto Prop::phaseInternalEnergy(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseInternalEnergy has not been implemented yet."); // TODO: Implement method Prop::phaseInternalEnergy.
}

auto Prop::phaseHelmholtzEnergy(const ChemicalSystem& system, const String& phase) -> Prop
{
    errorif(true, "Prop::phaseHelmholtzEnergy has not been implemented yet."); // TODO: Implement method Prop::phaseHelmholtzEnergy.
}

auto Prop::temperature(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::temperature has not been implemented yet."); // TODO: Implement method Prop::temperature.
}

auto Prop::pressure(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::pressure has not been implemented yet."); // TODO: Implement method Prop::pressure.
}

auto Prop::amount(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::amount has not been implemented yet."); // TODO: Implement method Prop::amount.
}

auto Prop::mass(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::mass has not been implemented yet."); // TODO: Implement method Prop::mass.
}

auto Prop::volume(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::volume has not been implemented yet."); // TODO: Implement method Prop::volume.
}

auto Prop::gibbsEnergy(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::gibbsEnergy has not been implemented yet."); // TODO: Implement method Prop::gibbsEnergy.
}

auto Prop::enthalpy(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::enthalpy has not been implemented yet."); // TODO: Implement method Prop::enthalpy.
}

auto Prop::entropy(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::entropy has not been implemented yet."); // TODO: Implement method Prop::entropy.
}

auto Prop::internalEnergy(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::internalEnergy has not been implemented yet."); // TODO: Implement method Prop::internalEnergy.
}

auto Prop::helmholtzEnergy(const ChemicalSystem& system) -> Prop
{
    errorif(true, "Prop::helmholtzEnergy has not been implemented yet."); // TODO: Implement method Prop::helmholtzEnergy.
}

} // namespace Reaktoro
