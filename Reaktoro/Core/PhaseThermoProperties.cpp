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
//#include "PhaseThermoProperties.hpp"
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
//struct PhaseThermoProperties::Impl
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
//    /// The results of the evaluation of the PhaseThermoModel function of the phase.
//    ThermoModelResult tres;
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
//};
//
//PhaseThermoProperties::PhaseThermoProperties()
//: pimpl(new Impl())
//{}
//
//PhaseThermoProperties::PhaseThermoProperties(const Phase& phase)
//: pimpl(new Impl(phase))
//{}
//
//PhaseThermoProperties::PhaseThermoProperties(const PhaseThermoProperties& other)
//: pimpl(new Impl(*other.pimpl))
//{}
//
//PhaseThermoProperties::~PhaseThermoProperties()
//{}
//
//auto PhaseThermoProperties::operator=(PhaseThermoProperties other) -> PhaseThermoProperties&
//{
//    pimpl = std::move(other.pimpl);
//    return *this;
//}
//
//auto PhaseThermoProperties::update(double T,double P) -> void
//{
//    pimpl->update(T, P);
//}
//
//auto PhaseThermoProperties::temperature() const -> double
//{
//    return pimpl->T.val;
//}
//
//auto PhaseThermoProperties::pressure() const -> double
//{
//    return pimpl->P.val;
//}
//
//auto PhaseThermoProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarGibbsEnergies();
//}
//
//auto PhaseThermoProperties::standardPartialMolarEnthalpies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarEnthalpies();
//}
//
//auto PhaseThermoProperties::standardPartialMolarVolumes() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarVolumes();
//}
//
//auto PhaseThermoProperties::standardPartialMolarEntropies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarEntropies();
//}
//
//auto PhaseThermoProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarInternalEnergies();
//}
//
//auto PhaseThermoProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHelmholtzEnergies();
//}
//
//auto PhaseThermoProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHeatCapacitiesConstP();
//}
//
//auto PhaseThermoProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
//{
//    return pimpl->standardPartialMolarHeatCapacitiesConstV();
//}
//
//} // namespace Reaktoro
