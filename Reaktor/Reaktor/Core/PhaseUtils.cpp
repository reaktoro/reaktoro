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
#include <Reaktor/Common/ScalarResult.hpp>
#include <Reaktor/Common/VectorResult.hpp>
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

auto chemicalPotentials(const Phase& phase, double T, double P) -> Vector
{
	const unsigned size = numSpecies(phase);
    Vector res(size);
    for(unsigned i = 0; i < size; ++i)
        res[i] = chemicalPotential(phase.species(i), T, P);
    return res;
}

auto molarFractions(const Phase& phase, const SubVector& n) -> Vector
{
	const unsigned size = n.n_rows;
	Assert(size == numSpecies(phase),
		"The dimension of the vector `n` is not the same "
		"as the number of species in the phase.");
    const double ntotal = arma::sum(n);
    if(ntotal == 0.0)
    	return zeros(size);
    return n/ntotal;
}

auto concentrations(const Phase& phase, const SubVector& n) -> Vector
{
    return phase.concentration()(n);
}

auto activity(const Phase& phase, const Index& i, double T, double P, const SubVector& n) -> ScalarResult
{
    return phase.activities()[i](T, P, n);
}

auto activities(const Phase& phase, double T, double P, const SubVector& n) -> VectorResult
{
	const unsigned size = numSpecies(phase);
	VectorResult res(size, size);
	for(unsigned i = 0; i < size; ++i)
		res.row(i) = activity(phase, i, T, P, n);
    return res;
}

auto names(const std::vector<Phase>& phases) -> std::vector<std::string>
{
	std::vector<std::string> names;
	names.reserve(phases.size());
	for(const Phase& phase : phases)
		names.push_back(phase.name());
	return names;
}

} /* namespace Reaktor */
