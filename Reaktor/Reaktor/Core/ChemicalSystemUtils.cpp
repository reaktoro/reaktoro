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

#include "ChemicalSystemUtils.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/ScalarResult.hpp>
#include <Reaktor/Common/VectorResult.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/PhaseUtils.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/SpeciesUtils.hpp>

namespace Reaktor {

auto numPhases(const ChemicalSystem& system) -> unsigned
{
	return system.phases().size();
}

auto numSpecies(const ChemicalSystem& system) -> unsigned
{
	return system.species().size();
}

auto numElements(const ChemicalSystem& system) -> unsigned
{
	return system.elements().size();
}

auto containsSpecies(const ChemicalSystem& system, const std::string& species) -> bool
{
	return indexSpecies(system, species) < numSpecies(system);
}

auto containsElement(const ChemicalSystem& system, const std::string& element) -> bool
{
	return indexElement(system, element) < numElements(system);
}

auto containsPhase(const ChemicalSystem& system, const std::string& phase) -> bool
{
	return indexPhase(system, phase) < numPhases(system);
}

auto indexElement(const ChemicalSystem& system, const std::string& name) -> Index
{
	const auto& begin = system.elements().begin();
	const auto& end = system.elements().end();
	return std::find(begin, end, name) - begin;
}

auto indicesElements(const ChemicalSystem& system, const std::vector<std::string>& names) -> Indices
{
	Indices indices;
	indices.reserve(names.size());
	for(const std::string& name : names)
		indices.push_back(indexElement(system, name));
	return indices;
}

auto indexSpecies(const ChemicalSystem& system, const std::string& name) -> Index
{
	const auto& begin = system.species().begin();
	const auto& end = system.species().end();
	const auto same_name = [&](const Species& species) { return species.name() == name; };
	return std::find_if(begin, end, same_name) - begin;
}

auto indicesSpecies(const ChemicalSystem& system, const std::vector<std::string>& names) -> Indices
{
	Indices indices;
	indices.reserve(names.size());
	for(const std::string& name : names)
		indices.push_back(indexSpecies(system, name));
	return indices;
}

auto indexPhase(const ChemicalSystem& system, const std::string& name) -> Index
{
	const auto& begin = system.phases().begin();
	const auto& end = system.phases().end();
	const auto same_name = [&](const Phase& phase) { return phase.name() == name; };
	return std::find_if(begin, end, same_name) - begin;
}

auto indicesPhases(const ChemicalSystem& system, const std::vector<std::string>& names) -> Indices
{
	Indices indices;
	indices.reserve(names.size());
	for(const std::string& name : names)
		indices.push_back(indexPhase(system, name));
	return indices;
}

auto indexFirstSpeciesInPhase(const ChemicalSystem& system, const Index& iphase) -> Index
{
	if(iphase < numPhases(system))
	{
		unsigned offset = 0;
		for(unsigned i = 0; i < iphase; ++i)
			offset += numSpecies(system.phase(i));
		return offset;
	}
	else return numSpecies(system);
}

auto indexLastSpeciesInPhase(const ChemicalSystem& system, const Index& iphase) -> Index
{
	return indexFirstSpeciesInPhase(system, iphase) +
		numSpecies(system.phase(iphase));
}

auto indicesElementsInSpecies(const ChemicalSystem& system, const Index& ispecies) -> Indices
{
	const Species& species = system.species(ispecies);
	return indicesElements(system, species.elements());
}

auto indicesSpeciesInPhase(const ChemicalSystem& system, const Index& iphase) -> Indices
{
	if(iphase < numPhases(system))
	{
		const unsigned num_species = numSpecies(system.phase(iphase));
		Indices indices(num_species);
		const Index first = indexFirstSpeciesInPhase(system, iphase);
		std::iota(indices.begin(), indices.end(), first);
		return indices;
	}
	else return Indices();
}

auto indicesSpeciesWithElement(const ChemicalSystem& system, const Index& ielement) -> Indices
{
	if(ielement < numElements(system))
	{
		std::set<Index> indices;
		for(unsigned i = 0; i < numSpecies(system); ++i)
			if(containsElement(system.species(i), system.element(ielement)))
				indices.insert(i);
		return Indices(indices.begin(), indices.end());
	}
	else return Indices();
}

auto indexPhaseWithSpecies(const ChemicalSystem& system, const Index& ispecies) -> Index
{
	if(ispecies < numSpecies(system))
	{
		for(unsigned iphase = 0; iphase < numPhases(system); ++iphase)
			if(ispecies < indexLastSpeciesInPhase(system, iphase))
				return iphase;
	}
	return numPhases(system);
}

auto localIndexSpecies(const ChemicalSystem& system, const Index& ispecies) -> Index
{
	const Index iphase = indexPhaseWithSpecies(system, ispecies);
	const Index ifirst = indexFirstSpeciesInPhase(system, iphase);
	return ispecies - ifirst;
}

auto chemicalPotentials(const ChemicalSystem& system, double T, double P) -> Vector
{
	Vector res(numSpecies(system));
	Index ifirst = 0;
	for(const Phase& phase : system.phases())
	{
		const unsigned size = numSpecies(phase);
		res.subvec(ifirst, ifirst + size) = chemicalPotentials(phase, T, P);
		ifirst += size;
	}
	return res;
}

auto molarFractions(const ChemicalSystem& system, const Vector& n) -> Vector
{
	Vector res(numSpecies(system));
	Index ifirst = 0;
	for(const Phase& phase : system.phases())
	{
		const Index size = numSpecies(phase);
		const Index ilast = ifirst + size;
		const SubVector nphase = n.subvec(ifirst, ilast);
		res.subvec(ifirst, ilast) = molarFractions(phase, nphase);
		ifirst += size;
	}
	return res;
}

auto concentrations(const ChemicalSystem& system, const Vector& n) -> Vector
{
	Vector res(numSpecies(system));
	Index ifirst = 0;
	for(const Phase& phase : system.phases())
	{
		const Index size = numSpecies(phase);
		const Index ilast = ifirst + size;
		const SubVector nphase = n.subvec(ifirst, ilast);
		res.subvec(ifirst, ilast) = concentrations(phase, nphase);
		ifirst += size;
	}
	return res;
}

auto activity(const ChemicalSystem& system, const Index& ispecies, double T, double P, const Vector& n) -> ScalarResult
{
	const Index iphase = indexPhaseWithSpecies(system, ispecies);
	const Index ilocalspecies = localIndexSpecies(system, ispecies);
	const Phase& phase = system.phase(iphase);
	const SubVector nphase = subvector(system, iphase, n);
	return activity(phase, ilocalspecies, T, P, nphase);
}

auto activities(const ChemicalSystem& system, double T, double P, const Vector& n) -> VectorResult
{
	VectorResult res(numSpecies(system));
	Index ifirst = 0;
	for(const Phase& phase : system.phases())
	{
		const Index size = numSpecies(phase);
		const Index ilast = ifirst + size;
		const SubVector nphase = n.subvec(ifirst, ilast);
		for(unsigned i = 0; i < size; ++i)
			res.row(i) = activity(phase, i, T, P, nphase);
		ifirst += size;
	}
	return res;
}

auto formulaMatrix(const ChemicalSystem& system) -> Matrix
{
	const auto& elements = system.elements();
	const auto& species = system.species();
	const auto& num_elements = elements.size();
	const auto& num_species = species.size();
	Matrix res(num_elements, num_species);
	for(unsigned i = 0; i < num_elements; ++i)
		for(unsigned j = 0; j < num_species; ++j)
			res(i, j) = elementAtoms(species[j], elements[i]);
	return res;
}

auto subvector(const ChemicalSystem& system, const Index& iphase, const Vector& vec) -> SubVector
{
	const Index first = indexFirstSpeciesInPhase(system, iphase);
	const Index last = indexLastSpeciesInPhase(system, iphase);
	return vec.subvec(first, last);
}

} /* namespace Reaktor */


