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

#include "ChemicalSystem.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {
namespace {

auto collectSpecies(const std::vector<Phase>& phases) -> std::vector<Species>
{
	std::vector<Species> species;
	for(const Phase& phase : phases)
		species.insert(species.end(),
			phase.species().begin(),
				phase.species().end());
	return species;
}

auto collectElements(const std::vector<Species>& species) -> std::vector<std::string>
{
	std::set<std::string> elements;
	for(const Species& iter : species)
		elements.insert(iter.elements().begin(), iter.elements().end());
	return std::vector<std::string>(elements.begin(), elements.end());
}

}

struct ChemicalSystem::Impl
{
    Impl()
    {}

    Impl(const std::vector<Phase>& phases)
    : phases(phases),
      species(collectSpecies(phases)),
      elements(collectElements(species))
    {}

    /// The phases in the system
    std::vector<Phase> phases;

    /// The chemical species in the system
    std::vector<Species> species;

    /// The chemical elements in the system
    std::vector<std::string> elements;
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

auto ChemicalSystem::phases() const -> const std::vector<Phase>&
{
	return pimpl->phases;
}

auto ChemicalSystem::phase(const Index& i) const -> const Phase&
{
	return pimpl->phases[i];
}

auto ChemicalSystem::species() const -> const std::vector<Species>&
{
	return pimpl->species;
}

auto ChemicalSystem::species(const Index& i) const -> const Species&
{
	return pimpl->species[i];
}

auto ChemicalSystem::elements() const -> const std::vector<std::string>&
{
    return pimpl->elements;
}

auto ChemicalSystem::element(const Index& i) const -> const std::string&
{
    return pimpl->elements[i];
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    out << "TODO" << std::endl; // todo implement operator<< for ChemicalSystem

	return out;
}

} /* namespace Reaktor */
