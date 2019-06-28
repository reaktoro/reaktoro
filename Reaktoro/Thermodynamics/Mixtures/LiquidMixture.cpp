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

#include "LiquidMixture.hpp"

namespace Reaktoro {

	LiquidMixture::LiquidMixture()
		: GeneralMixture<LiquidSpecies>()
	{}

	LiquidMixture::LiquidMixture(const std::vector<LiquidSpecies>& species)
		: GeneralMixture<LiquidSpecies>(species)
	{}

	LiquidMixture::~LiquidMixture()
	{}

	auto LiquidMixture::state(Temperature T, Pressure P, VectorConstRef n) const -> LiquidMixtureState
	{
		LiquidMixtureState res;
		res.T = T;
		res.P = P;
		res.x = moleFractions(n);
		return res;
	}

} // namespace Reaktoro
