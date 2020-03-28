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
//
//#pragma once
//
//// Reaktoro includes
////#include <Reaktoro/Math/Matrix.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class ChemicalState;
//class ChemicalSystem;
//class Partition;
//
///// A type that contains the values of a scalar field and its derivatives.
//class EquilibriumCompositionProblem
//{
//public:
//    /// Construct a default EquilibriumCompositionProblem instance.
//    EquilibriumCompositionProblem();
//
//    /// Construct a custom EquilibriumCompositionProblem instance.
//    EquilibriumCompositionProblem(const ChemicalSystem& system);
//
//    /// Construct a copy of a EquilibriumCompositionProblem instance.
//    EquilibriumCompositionProblem(const EquilibriumCompositionProblem& other);
//
//    /// Destroy this instance.
//    virtual ~EquilibriumCompositionProblem();
//
//    /// Construct a copy of a EquilibriumCompositionProblem instance.
//    auto operator=(EquilibriumCompositionProblem other) -> EquilibriumCompositionProblem&;
//
//    /// Return the chemical system.
//    auto system() const -> const ChemicalSystem&;
//
//    /// Return the partition of the chemical system.
//    auto partition() const -> const Partition&;
//
//    /// Set the partition of the chemical system.
//    auto setPartition(const Partition& partition) -> void;
//
//    /// Set the temperature for the equilibrium calculation.
//    /// @param value The temperature value.
//    /// @param units The temperature units.
//    auto setTemperature(double value, std::string units) -> void;
//
//    /// Set the pressure for the equilibrium calculation.
//    /// @param value The pressure value.
//    /// @param units The pressure units.
//    auto setPressure(double value, std::string units) -> void;
//
//    /// Set the composition of the aqueous phase using molalities of compounds.
//    /// The compounds and their molalities are separated by semicollon.
//    /// The following describes how to set the composition of an aqueous phase
//    /// with 1 molal of NaCl and 1 mmolal MgCl2:
//    /// ~~~
//    /// EquilibriumCompositionProblem composition(system);
//    /// composition.aqueous("1 molal NaCl; 1 mmolal MgCl2");
//    /// ~~~
//    auto setAqueousComposition(std::string molalities) -> void;
//
//    /// Set the composition of the gaseous phase using mole fractions of compounds.
//    /// The compounds and their mole fractions are separated by semicollon.
//    /// The following describes how to set the composition of a gas phase
//    /// with 70% N2, 20% O2, and 10% CO2 (molar percentage):
//    /// ~~~
//    /// EquilibriumCompositionProblem composition(system);
//    /// composition.gaseous("0.70 N2; 0.20 O2; 0.10 CO2");
//    /// ~~~
//    auto setGaseousComposition(std::string molarfractions) -> void;
//
//    /// Set the volume fractions of the solid phases.
//    /// The composition of the solid part of the system is defined using
//    /// volume fractions of each solid phase. The volume fraction of a solid
//    /// phase is defined as the volume of that phase divided by total solid volume.
//    /// The following describes how to set the volume fractions of solid
//    /// phases `Calcite` and `Quartz`.
//    /// ~~~
//    /// EquilibriumCompositionProblem composition(system);
//    /// composition.solid("0.10 Calcite; 0.90 Quartz");
//    /// ~~~
//    auto setSolidComposition(std::string volumefractions) -> void;
//
//    /// Set the saturation of the aqueous fluid.
//    /// The saturation of the aqueous fluid is defined as the ratio
//    /// of its volume and the total fluid volume.
//    auto setAqueousSaturation(double value) -> void;
//
//    /// Set the saturation of the gaseous fluid.
//    /// The saturation of the gaseous fluid is defined as the ratio
//    /// of its volume and the total fluid volume.
//    auto setGaseousSaturation(double value) -> void;
//
//    /// Set the porosity of the solid matrix.
//    /// The porosity is defined as the total fluid volume divided by total volume.
//    auto setPorosity(double value) -> void;
//
//    /// Convert this EquilibriumCompositionProblem instance into an EquilibriumProblem instance.
//    /// This conversion is needed to calculate the equilibrium state of both fluid and
//    /// solid phases using their given compositions and volume conditions.
//    /// Note that the calculated equilibrium state will satisfy the given fluid phase
//    /// saturations and solid matrix porosity. The internal equilibrium composition of
//    /// each phase might differ from those provided.
//    /// For example, assume the aqueous and gaseous phases are set as:
//    /// ~~~
//    /// composition.setAqueousFluid("1 molal NaCl");
//    /// composition.setGaseousFluid("0.95 CO2; 0.05 O2");
//    /// ChemicalState state = equilibrate(composition);
//    /// ~~~
//    /// When both phases are equilibrated, enouth amount of gas with prescribed
//    /// composition will be added in the system to satisfy the saturation of the
//    /// gaseous phase. As a result, the aqueous phase will become saturated with
//    /// both CO2 and O2. Thus, its final composition will contain a saturated
//    /// molality of CO2 and O2 in addition to NaCl.
//    operator EquilibriumInverseProblem();
//
//private:
//    struct Impl;
//
//    std::unique_ptr<Impl> pimpl;
//};
//
//} // namespace Reaktoro
//
