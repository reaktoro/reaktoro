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

#include "MultiphaseUtils.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/PhaseUtils.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/SpeciesUtils.hpp>

namespace Reaktor {

auto numElements(const Multiphase& multiphase) -> unsigned
{
    return multiphase.elements().size();
}

auto numSpecies(const Multiphase& multiphase) -> unsigned
{
    return multiphase.species().size();
}

auto numPhases(const Multiphase& multiphase) -> unsigned
{
    return multiphase.phases().size();
}

auto containsElement(const Multiphase& multiphase, const std::string& element) -> bool
{
    return indexElement(multiphase, element) < numElements(multiphase);
}

auto containsSpecies(const Multiphase& multiphase, const std::string& species) -> bool
{
    return indexSpecies(multiphase, species) < numSpecies(multiphase);
}

auto containsPhase(const Multiphase& multiphase, const std::string& phase) -> bool
{
    return indexPhase(multiphase, phase) < numPhases(multiphase);
}

auto indexElement(const Multiphase& multiphase, const std::string& element) -> Index
{
    const auto& begin = multiphase.elements().begin();
    const auto& end = multiphase.elements().end();
    return std::find(begin, end, element) - begin;
}

auto indicesElements(const Multiphase& multiphase, const std::vector<std::string>& names) -> Indices
{
    Indices indices;
    indices.reserve(names.size());
    for(const std::string& name : names)
        indices.push_back(indexElement(multiphase, name));
    return indices;
}

auto indexSpecies(const Multiphase& multiphase, const std::string& name) -> Index
{
    const auto& begin = multiphase.species().begin();
    const auto& end = multiphase.species().end();
    const auto sameName = [&](const Species& species) { return species.name() == name; };
    return std::find_if(begin, end, sameName) - begin;
}

auto indicesSpecies(const Multiphase& multiphase, const std::vector<std::string>& names) -> Indices
{
    Indices indices;
    indices.reserve(names.size());
    for(const std::string& name : names)
        indices.push_back(indexSpecies(multiphase, name));
    return indices;
}

auto indexPhase(const Multiphase& multiphase, const std::string& name) -> Index
{
    const auto& begin = multiphase.phases().begin();
    const auto& end = multiphase.phases().end();
    const auto sameName = [&](const Phase& phase) { return phase.name() == name; };
    return std::find_if(begin, end, sameName) - begin;
}

auto indicesPhases(const Multiphase& multiphase, const std::vector<std::string>& phases) -> Indices
{
    Indices indices;
    indices.reserve(phases.size());
    for(const std::string& name : phases)
        indices.push_back(indexPhase(multiphase, name));
    return indices;
}

auto indexFirstSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Index
{
    if(iphase < numPhases(multiphase))
    {
        unsigned offset = 0;
        for(unsigned i = 0; i < iphase; ++i)
            offset += numSpecies(multiphase.phases()[i]);
        return offset;
    }
    else return numSpecies(multiphase);
}

auto indexLastSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Index
{
    return indexFirstSpeciesInPhase(multiphase, iphase) +
        numSpecies(multiphase.phases()[iphase]);
}

auto indicesElementsInSpecies(const Multiphase& multiphase, const Index& ispecies) -> Indices
{
    const Species& species = multiphase.species()[ispecies];
    return indicesElements(multiphase, elementNames(species));
}

auto indicesElementsInSpecies(const Multiphase& multiphase, const Indices& ispecies) -> Indices
{
	std::set<Index> ielements;
	for(const Index& i : ispecies)
	{
		Indices tmp = indicesElementsInSpecies(multiphase, i);
		ielements.insert(tmp.begin(), tmp.end());
	}
	return Indices(ielements.begin(), ielements.end());
}

auto indicesSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Indices
{
    if(iphase < numPhases(multiphase))
    {
        const unsigned num_species = numSpecies(multiphase.phases()[iphase]);
        Indices indices(num_species);
        const Index first = indexFirstSpeciesInPhase(multiphase, iphase);
        std::iota(indices.begin(), indices.end(), first);
        return indices;
    }
    else return Indices();
}

auto indicesSpeciesWithElement(const Multiphase& multiphase, const Index& ielement) -> Indices
{
    if(ielement < numElements(multiphase))
    {
        std::set<Index> indices;
        for(unsigned i = 0; i < numSpecies(multiphase); ++i)
            if(containsElement(multiphase.species()[i], multiphase.elements()[ielement]))
                indices.insert(i);
        return Indices(indices.begin(), indices.end());
    }
    else return Indices();
}

auto indexPhaseWithSpecies(const Multiphase& multiphase, const Index& ispecies) -> Index
{
    if(ispecies < numSpecies(multiphase))
    {
        for(unsigned iphase = 0; iphase < numPhases(multiphase); ++iphase)
            if(ispecies < indexLastSpeciesInPhase(multiphase, iphase))
                return iphase;
    }
    return numPhases(multiphase);
}

auto indicesPhasesWithSpecies(const Multiphase& multiphase, const Indices& ispecies) -> Indices
{
	std::set<Index> iphases;
	for(const Index& i : ispecies)
		iphases.insert(indexPhaseWithSpecies(multiphase, i));
	return Indices(iphases.begin(), iphases.end());
}

auto localIndexSpecies(const Multiphase& multiphase, const Index& ispecies) -> Index
{
    const Index iphase = indexPhaseWithSpecies(multiphase, ispecies);
    const Index ifirst = indexFirstSpeciesInPhase(multiphase, iphase);
    return ispecies - ifirst;
}

auto indexMapSpeciesToElements(const Multiphase& multiphase) -> std::vector<Indices>
{
	const unsigned num_species = numSpecies(multiphase);
	std::vector<Indices> map(num_species);
	for(unsigned i = 0; i < num_species; ++i)
		map[i] = indicesElementsInSpecies(multiphase, i);
	return map;
}

auto indexMapElementToSpecies(const Multiphase& multiphase) -> std::vector<Indices>
{
	const unsigned num_elements = numElements(multiphase);
	std::vector<Indices> map(num_elements);
	for(unsigned i = 0; i < num_elements; ++i)
		map[i] = indicesSpeciesWithElement(multiphase, i);
	return map;
}

auto indexMapPhaseToSpecies(const Multiphase& multiphase) -> std::vector<Indices>
{
	const unsigned num_phases = numPhases(multiphase);
	std::vector<Indices> map(num_phases);
	for(unsigned i = 0; i < num_phases; ++i)
		map[i] = indicesSpeciesInPhase(multiphase, i);
	return map;
}

auto indexMapSpeciesToPhase(const Multiphase& multiphase) -> Indices
{
	const unsigned num_species = numSpecies(multiphase);
	Indices map(num_species);
	for(unsigned i = 0; i < num_species; ++i)
		map[i] = indexPhaseWithSpecies(multiphase, i);
	return map;
}

auto formulaMatrix(const Multiphase& multiphase) -> Matrix
{
    const auto& elements = multiphase.elements();
    const auto& species = multiphase.species();
    const auto& numElements = elements.size();
    const auto& numSpecies = species.size();
    Matrix res(numElements, numSpecies);
    for(unsigned i = 0; i < numElements; ++i)
        for(unsigned j = 0; j < numSpecies; ++j)
            res(i, j) = elementAtoms(species[j], elements[i]);
    return res;
}

auto subvector(const Multiphase& multiphase, const Index& iphase, const Vector& vec) -> VectorView
{
    const Index first = indexFirstSpeciesInPhase(multiphase, iphase);
    const Index last = indexLastSpeciesInPhase(multiphase, iphase);
    return vec.subvec(first, last);
}

auto volumes(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return volumes(multiphase.species(), T, P);
}

auto entropies(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return entropies(multiphase.species(), T, P);
}

auto helmholtzEnergies(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return helmholtzEnergies(multiphase.species(), T, P);
}

auto internalEnergies(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return internalEnergies(multiphase.species(), T, P);
}

auto enthalpies(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return enthalpies(multiphase.species(), T, P);
}

auto gibbsEnergies(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return gibbsEnergies(multiphase.species(), T, P);
}

auto heatCapacitiesCp(const Multiphase& multiphase, double T, double P) -> ThermoVector
{
    return heatCapacitiesCp(multiphase.species(), T, P);
}

auto molarFractions(const Multiphase& multiphase, const Vector& n) -> ThermoVector
{
    const unsigned nspecies = numSpecies(multiphase);
    ThermoVector res(nspecies, nspecies);
    Index ibegin = 0;
    for(const Phase& phase : multiphase.phases())
    {
        const Index size = numSpecies(phase);
        const Index iend = ibegin + size;
        const Vector nphase = n.subvec(ibegin, iend);
        for(unsigned i = 0; i < size; ++i)
            res.rows(ibegin, iend) = molarFractions(phase, nphase);
        ibegin += size;
    }
    return res;
}

auto concentrations(const Multiphase& multiphase, const Vector& n) -> ThermoVector
{
    const unsigned nspecies = numSpecies(multiphase);
    ThermoVector res(nspecies, nspecies);
    Index ibegin = 0;
    for(const Phase& phase : multiphase.phases())
    {
        const Index size = numSpecies(phase);
        const Index iend = ibegin + size;
        const Vector nphase = n.subvec(ibegin, iend);
        for(unsigned i = 0; i < size; ++i)
            res.rows(ibegin, iend) = concentrations(phase, nphase);
        ibegin += size;
    }
    return res;
}

auto activities(const Multiphase& multiphase, double T, double P, const Vector& n) -> ThermoVector
{
    const unsigned nspecies = numSpecies(multiphase);
    ThermoVector res(nspecies, nspecies);
    Index ibegin = 0;
    for(const Phase& phase : multiphase.phases())
    {
        const Index size = numSpecies(phase);
        const Index iend = ibegin + size;
        const Vector nphase = n.subvec(ibegin, iend);
        for(unsigned i = 0; i < size; ++i)
            res.rows(ibegin, iend) = activities(phase, T, P, nphase);
        ibegin += size;
    }
    return res;
}

} // namespace Reaktor

