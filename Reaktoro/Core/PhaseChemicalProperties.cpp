//// Reaktoro is a unified framework for modeling chemically reactive phases.
////
//// Copyright (C) 2014-2015 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#include "PhaseChemicalProperties.hpp"
//
//// Reaktoro includes
//#include <Reaktoro/Common/ChemicalScalar.hpp>
//#include <Reaktoro/Common/Constants.hpp>
//#include <Reaktoro/Common/ThermoScalar.hpp>
//#include <Reaktoro/Core/Phase.hpp>
//#include <Reaktoro/Core/Utils.hpp>
//#include <Reaktoro/Thermodynamics/Models/ChemicalModel.hpp>
//#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>
//
//namespace Reaktoro {
//
//struct PhaseChemicalProperties::Impl
//{
//    /// The phase
//    Phase phase;
//
//    /// The temperature of the phase (in units of K)
//    Temperature T;
//
//    /// The pressure of the phase (in units of Pa)
//    Pressure P;
//
//    /// The amounts of the species in the phase (in units of mol).
//    Vector n;
//
//    /// The molar fractions of the species in the phase (in units of mol/mol).
//    ChemicalVector x;
//
//    /// The results of the evaluation of the PhaseThermoModel function of the phase.
//    PhaseThermoModelResultConst tres;
//
//    /// The results of the evaluation of the PhaseChemicalModel function of the phase.
//    PhaseChemicalModelResultConst cres;
//
//    /// Construct a default Impl instance
//    Impl()
//    {}
//
//    /// Construct a Impl instance with given Phase
//    Impl(const Phase& phase)
//    : phase(phase)
//    {}
//
//    /// Update the thermodynamic properties of the phase.
//    auto update(double T_, double P_) -> void
//    {
//        // Set temperature and pressure
//        T = T_;
//        P = P_;
//
//        // Calculate the thermodynamic properties of the phase
//        auto tp = tres.map(0, phase.numSpecies());
//        phase.thermoModel()(tp, T_, P_);
//    }
//
//    /// Update the chemical properties of the phase.
//    auto update(double T_, double P_, const Vector& n_) -> void
//    {
//        // Set temperature, pressure, composition, and molar fractions
//        T = T_;
//        P = P_;
//        n = n_;
//        x = Reaktoro::molarFractions(n);
//
//        // Calculate the thermodynamic and chemical properties of the phase
//        auto tp = tres.map(0, phase.numSpecies());
//        auto cp = cres.map(0, 0, phase.numSpecies());
//        phase.thermoModel()(tp, T, P);
//        phase.chemicalModel()(cp, T, P, n_);
//    }
//
//    /// Return the molar fractions of the species.
//    auto molarFractions() const -> ChemicalVector
//    {
//        return x;
//    }
//
//    /// Return the ln activity coefficients of the species.
//    auto lnActivityCoefficients() const -> ChemicalVector
//    {
//        return cres.ln_activity_coefficients;
//    }
//
//    /// Return the ln activity constants of the species.
//    auto lnActivityConstants() const -> ThermoVector
//    {
//        return cres.ln_activity_constants;
//    }
//
//    /// Return the ln activities of the species.
//    auto lnActivities() const -> ChemicalVector
//    {
//        return cres.ln_activities;
//    }
//
//    /// Return the chemical potentials of the species (in units of J/mol).
//    auto chemicalPotentials() const -> ChemicalVector
//    {
//        const auto& R = universalGasConstant;
//        const auto& G = tres.standard_partial_molar_gibbs_energies;
//        const auto& lna = cres.ln_activities;
//        return G + R*T*lna;
//    }
//
//    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
//    auto standardPartialMolarGibbsEnergies() const -> ThermoVector
//    {
//        return tres.standard_partial_molar_gibbs_energies;
//    }
//
//    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
//    auto standardPartialMolarEnthalpies() const -> ThermoVector
//    {
//        return tres.standard_partial_molar_enthalpies;
//    }
//
//    /// Return the standard partial molar volumes of the species (in units of m3/mol).
//    auto standardPartialMolarVolumes() const -> ThermoVector
//    {
//        return tres.standard_partial_molar_volumes;
//    }
//
//    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
//    auto standardPartialMolarEntropies() const -> ThermoVector
//    {
//        const auto& G = standardPartialMolarGibbsEnergies();
//        const auto& H = standardPartialMolarEnthalpies();
//        return (H - G)/T;
//    }
//
//    /// Return the standard partial molar internal energies of the species (in units of J/mol).
//    auto standardPartialMolarInternalEnergies() const -> ThermoVector
//    {
//        const auto& H = standardPartialMolarEnthalpies();
//        const auto& V = standardPartialMolarVolumes();
//        return H - P*V;
//    }
//
//    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
//    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector
//    {
//        const auto& G = standardPartialMolarGibbsEnergies();
//        const auto& V = standardPartialMolarVolumes();
//        return G - P*V;
//    }
//
//    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
//    {
//        return tres.standard_partial_molar_heat_capacities_cp;
//    }
//
//    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
//    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
//    {
//        return tres.standard_partial_molar_heat_capacities_cv;
//    }
//
//    /// Return the molar Gibbs energy of the phase (in units of J/mol).
//    auto molarGibbsEnergy() const -> ChemicalScalar
//    {
//        ChemicalScalar res = sum(x % tres.standard_partial_molar_gibbs_energies);
//        res += cres.residual_molar_gibbs_energy[0];
//        return res;
//    }
//
//    /// Return the molar enthalpy of the phase (in units of J/mol).
//    auto molarEnthalpy() const -> ChemicalScalar
//    {
//        ChemicalScalar res = sum(x % tres.standard_partial_molar_enthalpies);
//        res += cres.residual_molar_enthalpy[0];
//        return res;
//    }
//
//    /// Return the molar volumes of the phase (in units of m3/mol).
//    auto molarVolume() const -> ChemicalScalar
//    {
//        if(cres.molar_volume[0] > 0.0)
//            return cres.molar_volume[0];
//        return sum(x % tres.standard_partial_molar_volumes);
//    }
//
//    /// Return the molar entropy of the phase (in units of J/(mol*K)).
//    auto molarEntropy() const -> ChemicalScalar
//    {
//        const auto& G = molarGibbsEnergy();
//        const auto& H = molarEnthalpy();
//        return (H - G)/T;
//    }
//
//    /// Return the molar internal energy of the phase (in units of J/mol).
//    auto molarInternalEnergy() const -> ChemicalScalar
//    {
//        const auto& H = molarEnthalpy();
//        const auto& V = molarVolume();
//        return H - P*V;
//    }
//
//    /// Return the molar Helmholtz energy of the phase (in units of J/mol).
//    auto molarHelmholtzEnergy() const -> ChemicalScalar
//    {
//        const auto& G = molarGibbsEnergy();
//        const auto& V = molarVolume();
//        return G - P*V;
//    }
//
//    /// Return the molar isobaric heat capacities of the phase (in units of J/(mol*K)).
//    auto molarHeatCapacityConstP() const -> ChemicalScalar
//    {
//        ChemicalScalar res = sum(x % tres.standard_partial_molar_heat_capacities_cp);
//        res += cres.residual_molar_heat_capacity_cp[0];
//        return res;
//    }
//
//    /// Return the molar isochoric heat capacities of the phase (in units of J/(mol*K)).
//    auto molarHeatCapacityConstV() const -> ChemicalScalar
//    {
//        ChemicalScalar res = sum(x % tres.standard_partial_molar_heat_capacities_cv);
//        res += cres.residual_molar_heat_capacity_cv[0];
//        return res;
//    }
//
//    /// Return the specific Gibbs energy of the phase (in units of J/kg).
//    auto specificGibbsEnergy() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarGibbsEnergy();
//    }
//
//    /// Return the specific enthalpy of the phase (in units of J/kg).
//    auto specificEnthalpy() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarEnthalpy();
//    }
//
//    /// Return the specific volumes of the phase (in units of m3/kg).
//    auto specificVolume() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarVolume();
//    }
//
//    /// Return the specific entropy of the phase (in units of J/(kg*K)).
//    auto specificEntropy() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarEntropy();
//    }
//
//    /// Return the specific internal energy of the phase (in units of J/kg).
//    auto specificInternalEnergy() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarInternalEnergy();
//    }
//
//    /// Return the specific Helmholtz energy of the phase (in units of J/kg).
//    auto specificHelmholtzEnergy() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarHelmholtzEnergy();
//    }
//
//    /// Return the specific isobaric heat capacities of the phase (in units of J/(kg*K)).
//    auto specificHeatCapacityConstP() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarHeatCapacityConstP();
//    }
//
//    /// Return the specific isochoric heat capacities of the phase (in units of J/(kg*K)).
//    auto specificHeatCapacityConstV() const -> ChemicalScalar
//    {
//        return amount()/mass() * molarHeatCapacityConstV();
//    }
//
//    /// Return the density of the phase (in units of kg/m3).
//    auto density() const -> ChemicalScalar
//    {
//        return mass()/volume();
//    }
//
//    /// Return the mass of the phase (in units of kg).
//    auto mass() const -> ChemicalScalar
//    {
//        auto nc = Reaktoro::composition(n);
//        auto mm = Reaktoro::molarMasses(phase.species());
//        return sum(mm % nc);
//    }
//
//    /// Return the amount of the phase (in units of mol).
//    auto amount() const -> ChemicalScalar
//    {
//        auto nc = Reaktoro::composition(n);
//        return sum(nc);
//    }
//
//    /// Return the volume of the phase (in units of m3).
//    auto volume() const -> ChemicalScalar
//    {
//        return amount() * molarVolume();
//    }
//};
//
//PhaseChemicalProperties::PhaseChemicalProperties()
//: pimpl(new Impl())
//{}
//
//PhaseChemicalProperties::PhaseChemicalProperties(const Phase& phase)
//: pimpl(new Impl(phase))
//{}
//
//PhaseChemicalProperties::PhaseChemicalProperties(const PhaseChemicalProperties& other)
//: pimpl(new Impl(*other.pimpl))
//{}
//
//PhaseChemicalProperties::~PhaseChemicalProperties()
//{}
//
//auto PhaseChemicalProperties::operator=(PhaseChemicalProperties other) -> PhaseChemicalProperties&
//{
//    pimpl = std::move(other.pimpl);
//    return *this;
//}
//
//auto PhaseChemicalProperties::update(double T,double P) -> void
//{
//    pimpl->update(T, P);
//}
//
//auto PhaseChemicalProperties::update(double T, double P, const Vector& n) -> void
//{
//    pimpl->update(T, P, n);
//}
//
//auto PhaseChemicalProperties::temperature() const -> double
//{
//    return pimpl->T.val;
//}
//
//auto PhaseChemicalProperties::pressure() const -> double
//{
//    return pimpl->P.val;
//}
//
//auto PhaseChemicalProperties::composition()const -> const Vector&
//{
//    return pimpl->n;
//}
//
//auto PhaseChemicalProperties::molarFractions() const -> ChemicalVector
//{
//    return pimpl->molarFractions();
//}
//
//auto PhaseChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
//{
//    return pimpl->lnActivityCoefficients();
//}
//
//auto PhaseChemicalProperties::lnActivityConstants() const -> ThermoVector
//{
//    return pimpl->lnActivityConstants();
//}
//
//auto PhaseChemicalProperties::lnActivities() const -> ChemicalVector
//{
//    return pimpl->lnActivities();
//}
//
//auto PhaseChemicalProperties::chemicalPotentials() const -> ChemicalVector
//{
//    return pimpl->chemicalPotentials();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarGibbsEnergies();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarEnthalpies();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarVolumes();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarEntropies();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarInternalEnergies();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHelmholtzEnergies();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHeatCapacitiesConstP();
//}
//
//auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHeatCapacitiesConstV();
//}
//
//auto PhaseChemicalProperties::molarGibbsEnergy() const -> ChemicalScalar
//{
//    return pimpl->molarGibbsEnergy();
//}
//
//auto PhaseChemicalProperties::molarEnthalpy() const -> ChemicalScalar
//{
//    return pimpl->molarEnthalpy();
//}
//
//auto PhaseChemicalProperties::molarVolume() const -> ChemicalScalar
//{
//    return pimpl->molarVolume();
//}
//
//auto PhaseChemicalProperties::molarEntropy() const -> ChemicalScalar
//{
//    return pimpl->molarEntropy();
//}
//
//auto PhaseChemicalProperties::molarInternalEnergy() const -> ChemicalScalar
//{
//    return pimpl->molarInternalEnergy();
//}
//
//auto PhaseChemicalProperties::molarHelmholtzEnergy() const -> ChemicalScalar
//{
//    return pimpl->molarHelmholtzEnergy();
//}
//
//auto PhaseChemicalProperties::molarHeatCapacityConstP() const -> ChemicalScalar
//{
//    return pimpl->molarHeatCapacityConstP();
//}
//
//auto PhaseChemicalProperties::molarHeatCapacityConstV() const -> ChemicalScalar
//{
//    return pimpl->molarHeatCapacityConstV();
//}
//
//auto PhaseChemicalProperties::specificGibbsEnergy() const -> ChemicalScalar
//{
//    return pimpl->specificGibbsEnergy();
//}
//
//auto PhaseChemicalProperties::specificEnthalpy() const -> ChemicalScalar
//{
//    return pimpl->specificEnthalpy();
//}
//
//auto PhaseChemicalProperties::specificVolume() const -> ChemicalScalar
//{
//    return pimpl->specificVolume();
//}
//
//auto PhaseChemicalProperties::specificEntropy() const -> ChemicalScalar
//{
//    return pimpl->specificEntropy();
//}
//
//auto PhaseChemicalProperties::specificInternalEnergy() const -> ChemicalScalar
//{
//    return pimpl->specificInternalEnergy();
//}
//
//auto PhaseChemicalProperties::specificHelmholtzEnergy() const -> ChemicalScalar
//{
//    return pimpl->specificHelmholtzEnergy();
//}
//
//auto PhaseChemicalProperties::specificHeatCapacityConstP() const -> ChemicalScalar
//{
//    return pimpl->specificHeatCapacityConstP();
//}
//
//auto PhaseChemicalProperties::specificHeatCapacityConstV() const -> ChemicalScalar
//{
//    return pimpl->specificHeatCapacityConstV();
//}
//
//auto PhaseChemicalProperties::density() const -> ChemicalScalar
//{
//    return pimpl->density();
//}
//
//auto PhaseChemicalProperties::mass() const -> ChemicalScalar
//{
//    return pimpl->mass();
//}
//
//auto PhaseChemicalProperties::amount() const -> ChemicalScalar
//{
//    return pimpl->amount();
//}
//
//auto PhaseChemicalProperties::volume() const -> ChemicalScalar
//{
//    return pimpl->volume();
//}
//
//} // namespace Reaktoro
