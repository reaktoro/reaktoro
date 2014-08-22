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

#include "Species.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>

namespace Reaktor {

SpeciesParams::SpeciesParams()
: molarMass(INFINITY), charge(INFINITY)
{}

struct Species::Impl
{
	Impl()
	{}

	Impl(const SpeciesParams& params) : params(params)
	{
		const auto initialisedName = not params.name.empty();
		const auto initialisedFormula = not params.formula.empty();
		const auto initialisedElements = not params.elements.empty();
		const auto initialisedCoefficients = not params.coefficients.empty();
		const auto initialisedChemicalPotential = bool(params.chemicalPotential);
		const auto initialisedActivity = bool(params.activity);
		const auto initialisedCharge = params.charge != INFINITY;
		const auto initialisedMolarMass = params.molarMass != INFINITY;
		const auto nonnegativeMolarMass = params.molarMass > 0;
		const auto numElements = params.elements.size();
		const auto numCoefficients = params.coefficients.size();

		Assert(initialisedName,
			"The `name` parameter of the species has not been initialised.");
		Assert(initialisedFormula,
			"The `formula` parameter of the species has not been initialised.");
		Assert(initialisedElements,
			"The `elements` parameter of the species has not been initialised.");
		Assert(initialisedCoefficients,
			"The `coefficients` parameter of the species has not been initialised.");
		Assert(initialisedChemicalPotential,
			"The `chemicalPotential` parameter of the species has not been initialised.");
		Assert(initialisedActivity,
			"The `activity` parameter of the species has not been initialised.");
		Assert(initialisedCharge,
			"The `charge` parameter of the species has not been initialised.");
		Assert(initialisedMolarMass,
			"The `molarMass` parameter of the species has not been initialised.");
		Assert(numElements == numCoefficients,
			"The number of elements is not the same as the number of coefficients.");
		Assert(nonnegativeMolarMass,
			"The 'molarMass' property of a species must be non-negative.");
	}

	SpeciesParams params;
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const SpeciesParams& params)
: pimpl(new Impl(params))
{}

auto Species::name() const -> const std::string&
{
	return pimpl->params.name;
}

auto Species::formula() const -> const std::string&
{
	return pimpl->params.formula;
}

auto Species::elements() const -> const std::vector<std::string>&
{
	return pimpl->params.elements;
}

auto Species::coefficients() const -> const std::vector<double>&
{
	return pimpl->params.coefficients;
}

auto Species::molarMass() const -> double
{
	return pimpl->params.molarMass;
}

auto Species::charge() const -> double
{
	return pimpl->params.charge;
}

auto Species::chemicalPotential() const -> const ChemicalPotential&
{
	return pimpl->params.chemicalPotential;
}

auto Species::activity() const -> const Activity&
{
	return pimpl->params.activity;
}

auto operator<<(std::ostream& out, const Species& species) -> std::ostream&
{
	// todo improve operator<< of Species
    out << "Species" << std::endl;
    out << "    name:         " << species.name() << std::endl;
    out << "    formula:      " << species.formula() << std::endl;
//    out << "    elements:     " << species.elements() << std::endl;
//    out << "    coefficients: " << species.coefficients() << std::endl;
    out << "    charge:       " << species.charge() << std::endl;
    out << "    molar-mass:   " << species.molarMass() << " mol/kg" << std::endl;

    return out;
}

} /* namespace Reaktor */
