// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PhaseUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/ThermoProperties.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/SpeciesUtils.hpp>

namespace Reaktor {

auto numSpecies(const Phase& phase) -> unsigned
{
    return phase.species().size();
}

auto indexSpecies(const Phase& phase, const std::string& name) -> Index
{
    const auto compare = [&](const Species& s) { return s.name() == name; };
    const auto& species = phase.species();
    return std::find_if(species.begin(), species.end(), compare) - species.begin();
}

auto containsSpecies(const Phase& phase, const std::string& name) -> bool
{
	return indexSpecies(phase, name) < numSpecies(phase);
}

auto names(const std::vector<Phase>& phases) -> std::vector<std::string>
{
	std::vector<std::string> names;
	names.reserve(phases.size());
	for(const Phase& phase : phases)
		names.push_back(phase.name());
	return names;
}

auto volumes(const Phase& phase, double T, double P) -> ThermoProperties
{
	return volumes(phase.species(), T, P);
}

auto entropies(const Phase& phase, double T, double P) -> ThermoProperties
{
	return entropies(phase.species(), T, P);
}

auto helmholtzEnergies(const Phase& phase, double T, double P) -> ThermoProperties
{
	return helmholtzEnergies(phase.species(), T, P);
}

auto internalEnergies(const Phase& phase, double T, double P) -> ThermoProperties
{
	return internalEnergies(phase.species(), T, P);
}

auto enthalpies(const Phase& phase, double T, double P) -> ThermoProperties
{
	return enthalpies(phase.species(), T, P);
}

auto gibbsEnergies(const Phase& phase, double T, double P) -> ThermoProperties
{
	return gibbsEnergies(phase.species(), T, P);
}

auto heatCapacitiesCp(const Phase& phase, double T, double P) -> ThermoProperties
{
	return heatCapacitiesCp(phase.species(), T, P);
}

auto molarFractions(const Phase& phase, const Vector& n) -> ThermoVector
{
	const unsigned nspecies = n.size();
	ThermoVector x = ThermoVector::zero(nspecies, nspecies);
	Assert(nspecies == numSpecies(phase),
		"The dimension of the vector `n` is not the same "
		"as the number of species in the phase.");
    const double ntotal = arma::sum(n);
    if(ntotal == 0.0)
    	return x;
    x.val = n/ntotal;
    for(unsigned i = 0; i < nspecies; ++i)
	{
		x.ddn(i, i) = 1.0/ntotal;
		x.ddn.row(i) -= x.val.row(i)/ntotal;
	}
    return x;
}

auto concentrations(const Phase& phase, const Vector& n) -> ThermoVector
{
    return phase.thermoModel().concentration(n);
}

auto activities(const Phase& phase, double T, double P, const Vector& n) -> ThermoVector
{
    return phase.thermoModel().activity(T, P, n);
}

} // namespace Reaktor
