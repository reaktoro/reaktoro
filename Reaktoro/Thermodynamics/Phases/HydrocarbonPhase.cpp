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

#include "HydrocarbonPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/HydrocarbonMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/LiquidChemicalModelCubicEOS.hpp>

namespace Reaktoro {

	struct HydrocarbonPhase::Impl
	{
		/// The hydrocarbon mixture instance
		HydrocarbonMixture mixture;

		/// Construct a default Impl instance
		Impl()
		{}

		/// Construct a custom Impl instance
		Impl(const HydrocarbonMixture& mixture)
			: mixture(mixture)
		{}
	};

	HydrocarbonPhase::HydrocarbonPhase()
		: Phase(), pimpl(new Impl())
	{}

	HydrocarbonPhase::HydrocarbonPhase(const HydrocarbonMixture& mixture)
		: pimpl(new Impl(mixture))
	{
		// Convert the HydrocarbonSpecies instances to Species instances
		std::vector<Species> species;
		for (const HydrocarbonSpecies& x : mixture.species())
			species.push_back(x);

		// Set the Phase attributes
		setName("Liquid");
		setType(PhaseType::Liquid);
		setSpecies(species);
		setChemicalModelPengRobinson();
	}

	auto HydrocarbonPhase::setChemicalModelVanDerWaals() -> HydrocarbonPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelVanDerWaals(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto HydrocarbonPhase::setChemicalModelRedlichKwong() -> HydrocarbonPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelRedlichKwong(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto HydrocarbonPhase::setChemicalModelSoaveRedlichKwong() -> HydrocarbonPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelSoaveRedlichKwong(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto HydrocarbonPhase::setChemicalModelPengRobinson() -> HydrocarbonPhase&
	{
		PhaseChemicalModel model = liquidChemicalModelPengRobinson(mixture());
		setChemicalModel(model);
		return *this;
	}

	auto HydrocarbonPhase::mixture() const -> const HydrocarbonMixture&
	{
		return pimpl->mixture;
	}

} // namespace Reaktoro
