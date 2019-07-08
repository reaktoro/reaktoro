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

#include "FluidPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/FluidChemicalModelCubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Models/FluidChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/FluidChemicalModelSpycherReed.hpp>
#include <Reaktoro/Thermodynamics/Models/FluidChemicalModelSpycherPruessEnnis.hpp>

namespace Reaktoro {

    struct FluidPhase::Impl
    {
        /// The fluid mixture instance
        FluidMixture mixture;

        /// Construct a default Impl instance
        Impl()
        {}

        /// Construct a custom Impl instance
        Impl(const FluidMixture& mixture)
            : mixture(mixture)
        {}

    };

    FluidPhase::FluidPhase()
        : Phase("Fluid", PhaseType::Fluid), pimpl(new Impl())
    {}

    FluidPhase::FluidPhase(const FluidMixture& mixture, std::string name, PhaseType type)
        : Phase(name, type), pimpl(new Impl(mixture))
    {
        // Convert the FluidSpecies instances to Species instances
        std::vector<Species> species;
        for (const FluidSpecies& x : mixture.species())
            species.push_back(x);

        // Set the Phase attributes
        setSpecies(species);
        setChemicalModelPengRobinson();
    }

    auto FluidPhase::setChemicalModelIdeal() -> FluidPhase&
    {
        Assert(type() == PhaseType::Gas, "Try to set Ideal model in a phase that is not defined as Gas, change PhaseType to Gas", "Ideal model can only be used to model a Gas phases");
        PhaseChemicalModel model = fluidChemicalModelIdeal(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelVanDerWaals() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelVanDerWaals(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelRedlichKwong() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelRedlichKwong(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelSoaveRedlichKwong() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelSoaveRedlichKwong(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelPengRobinson() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelPengRobinson(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelSpycherPruessEnnis() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelSpycherPruessEnnis(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::setChemicalModelSpycherReed() -> FluidPhase&
    {
        PhaseChemicalModel model = fluidChemicalModelSpycherReed(mixture());
        setChemicalModel(model);
        return *this;
    }

    auto FluidPhase::mixture() const -> const FluidMixture&
    {
        return pimpl->mixture;
    }

} // namespace Reaktoro
