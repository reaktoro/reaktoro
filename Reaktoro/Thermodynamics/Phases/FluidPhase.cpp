// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "FluidPhase.hpp"

// // Reaktoro includes
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Common/Index.hpp>
// #include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
// #include <Reaktoro/Thermodynamics/Models/FluidChemicalModelCubicEOS.hpp>
// #include <Reaktoro/Thermodynamics/Models/FluidChemicalModelIdeal.hpp>
// #include <Reaktoro/Thermodynamics/Models/FluidChemicalModelSpycherReed.hpp>
// #include <Reaktoro/Thermodynamics/Models/FluidChemicalModelSpycherPruessEnnis.hpp>

// namespace Reaktoro {

// struct FluidPhase::Impl
// {
//     /// The fluid mixture instance
//     GeneralMixture mixture;

//     /// Construct a default Impl instance
//     Impl()
//     {}

//     /// Construct a custom Impl instance
//     Impl(const GeneralMixture& mixture)
//         : mixture(mixture)
//     {}

// };

// FluidPhase::FluidPhase(const std::string& name, StateOfMatter type)
//     : Phase(name, type), pimpl(new Impl())
// {}

// FluidPhase::FluidPhase(const GeneralMixture& mixture, const std::string& name, StateOfMatter type)
//     : Phase(name, type), pimpl(new Impl(mixture))
// {
//     // Convert the FluidSpecies instances to Species instances
//     std::vector<Species> species;
//     for (const FluidSpecies& x : mixture.species())
//         species.push_back(x);

//     // Set the Phase attributes
//     setSpecies(species);
//     setChemicalModelPengRobinson();
// }

// auto FluidPhase::setChemicalModelIdeal() -> FluidPhase&
// {
//     Assert(type() == StateOfMatter::Gas, "Try to set Ideal model in a phase that is not defined as Gas, change StateOfMatter to Gas", "Ideal model can only be used to model a Gas phases");
//     ActivityPropsFn fn = fluidChemicalModelIdeal(mixture());
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelVanDerWaals() -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelCubicEOS(mixture(), type(), {CubicEOS::VanDerWaals});
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelRedlichKwong() -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelCubicEOS(mixture(), type(), {CubicEOS::RedlichKwong});
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelSoaveRedlichKwong() -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelCubicEOS(mixture(), type(), {CubicEOS::SoaveRedlichKwong});
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelPengRobinson(CubicEOS::Params params) -> FluidPhase&
// {
//     params.model = CubicEOS::PengRobinson;
//     ActivityPropsFn fn = fluidChemicalModelCubicEOS(mixture(), type(), params);
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelCubicEOS(CubicEOS::Params params) -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelCubicEOS(mixture(), type(), params);
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelSpycherPruessEnnis() -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelSpycherPruessEnnis(mixture());
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::setChemicalModelSpycherReed() -> FluidPhase&
// {
//     ActivityPropsFn fn = fluidChemicalModelSpycherReed(mixture());
//     setChemicalModel(fn);
//     return *this;
// }

// auto FluidPhase::mixture() const -> const GeneralMixture&
// {
//     return pimpl->mixture;
// }

// auto FluidPhase::mixture() -> GeneralMixture&
// {
//     return pimpl->mixture;
// }



// } // namespace Reaktoro
