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

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include <Reaktoro/Thermodynamics/Species/HydrocarbonSpecies.hpp>

namespace Reaktoro {

	/// A type used to describe the state of a hydrocarbon mixture
	struct HydrocarbonMixtureState : public MixtureState
	{};

	/// Provides a computational representation of a hydrocarbon mixture.
	/// The HydrocarbonMixture class is defined as a collection of HydrocarbonSpecies objects,
	/// representing, therefore, a mixture of hydrocarbon species. Its main purpose is to
	/// provide the necessary operations in the calculation of activities of hydrocarbon
	/// species.
	/// @see HydrocarbonSpecies
	/// @ingroup Mixtures
	class HydrocarbonMixture : public GeneralMixture<HydrocarbonSpecies>
	{
	public:
		/// Construct a default HydrocarbonMixture instance.
		HydrocarbonMixture();

		/// Construct a HydrocarbonMixture instance with given species.
		/// @param species The species that compose the hydrocarbon mixture
		explicit HydrocarbonMixture(const std::vector<HydrocarbonSpecies>& species);

		/// Destroy the HydrocarbonMixture instance.
		virtual ~HydrocarbonMixture();

		/// Calculate the state of the hydrocarbon mixture.
		/// @param T The temperature (in units of K)
		/// @param P The pressure (in units of Pa)
		/// @param n The molar amounts of the species in the mixture (in units of mol)
		auto state(Temperature T, Pressure P, VectorConstRef n) const->HydrocarbonMixtureState;
	};

} // namespace Reaktoro
