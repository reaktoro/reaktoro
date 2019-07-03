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

#include "LiquidPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/LiquidMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/LiquidChemicalModelCubicEOS.hpp>

namespace Reaktoro {

	struct LiquidPhase::Impl
	{
		/// The liquid mixture instance
		LiquidMixture mixture;

		/// Construct a default Impl instance
		Impl()
		{}

		/// Construct a custom Impl instance
		Impl(const LiquidMixture& mixture)
			: mixture(mixture)
		{}
	};

	LiquidPhase::LiquidPhase()
		: Phase(), pimpl(new Impl())
	{
        setName("Liquid");
        setType(PhaseType::Liquid);
    }

	LiquidPhase::LiquidPhase(const LiquidMixture& mixture)
		: pimpl(new Impl(mixture))
	{
		// Convert the LiquidSpecies instances to Species instances
		std::vector<Species> species;
		for (const LiquidSpecies& x : mixture.species())
			species.push_back(x);

		// Set the Phase attributes
		setName("Liquid");
		setType(PhaseType::Liquid);
		setSpecies(species);
		setChemicalModelPengRobinson();
	}

	auto LiquidPhase::setChemicalModelVanDerWaals() -> LiquidPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelVanDerWaals(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto LiquidPhase::setChemicalModelRedlichKwong() -> LiquidPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelRedlichKwong(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto LiquidPhase::setChemicalModelSoaveRedlichKwong() -> LiquidPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelSoaveRedlichKwong(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto LiquidPhase::setChemicalModelPengRobinson() -> LiquidPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelPengRobinson(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto LiquidPhase::mixture() const -> const LiquidMixture&
	{
		return pimpl->mixture;
	}

} // namespace Reaktoro
