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
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {
namespace {

auto formulaMatrix(const std::vector<Species>& species, const std::vector<Element>& elements) -> Matrix
{
    const auto& num_elements = elements.size();
    const auto& num_species = species.size();
    Matrix W(num_elements, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned j = 0; j < num_elements; ++j)
            W(j, i) = atoms(elements[j], species[i]);
    return W;
}

} // namespace

struct ChemicalSystem::Impl
{
    /// The data used to construct the chemical system
    ChemicalSystemData data;

    /// The list of species in the chemical system
    std::vector<Species> species;

    /// The list of elements in the chemical system
    std::vector<Element> elements;

    /// The formula matrix of the chemical system
    Matrix formula_matrix;

    Impl()
    {}

    Impl(const ChemicalSystemData& data)
    : data(data), species(Reaktor::species(data.phases)), elements(Reaktor::elements(species))
    {
        formula_matrix = Reaktor::formulaMatrix(species, elements);
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const ChemicalSystemData& data)
: pimpl(new Impl(data))
{}

auto ChemicalSystem::numElements() const -> unsigned
{
    return elements().size();
}

auto ChemicalSystem::numSpecies() const -> unsigned
{
    return species().size();
}

auto ChemicalSystem::numPhases() const -> unsigned
{
    return phases().size();
}

auto ChemicalSystem::element(Index index) const -> const Element&
{
    return elements()[index];
}

auto ChemicalSystem::element(std::string name) const -> const Element&
{
    return element(indexElement(name));
}

auto ChemicalSystem::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto ChemicalSystem::species(Index index) const -> const Species&
{
    return species()[index];
}

auto ChemicalSystem::species(std::string name) const -> const Species&
{
    return species(indexSpecies(name));
}

auto ChemicalSystem::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto ChemicalSystem::phase(Index index) const -> const Phase&
{
    return phases()[index];
}

auto ChemicalSystem::phase(std::string name) const -> const Phase&
{
    return phase(indexPhase(name));
}

auto ChemicalSystem::phases() const -> const std::vector<Phase>&
{
    return pimpl->data.phases;
}

auto ChemicalSystem::formulaMatrix() const -> const Matrix&
{
    return pimpl->formula_matrix;
}

auto ChemicalSystem::indexElement(std::string name) const -> Index
{
    return index(name, elements());
}

auto ChemicalSystem::indexElementWithError(std::string name) const -> Index
{
    const Index index = indexElement(name);

    Assert(index < numElements(),
        "Cannot get the index of the element " + name + ".",
        "There is no element called " + name + " in the chemical system.");

    return index;
}

auto ChemicalSystem::indexSpecies(std::string name) const -> Index
{
    return index(name, species());
}

auto ChemicalSystem::indexSpeciesWithError(std::string name) const -> Index
{
    const Index index = indexSpecies(name);

    Assert(index < numSpecies(),
        "Cannot get the index of the species " + name + ".",
        "There is no element called " + name + " in the chemical system.");

    return index;
}

auto ChemicalSystem::indexPhase(std::string name) const -> Index
{
    return index(name, phases());
}

auto ChemicalSystem::indexPhaseWithError(std::string name) const -> Index
{
    const Index index = indexPhase(name);

    Assert(index < numPhases(),
        "Cannot get the index of the phase " + name + ".",
        "There is no phase called " + name + " in the chemical system.");

    return index;
}

auto ChemicalSystem::indexPhaseWithSpecies(Index index) const -> Index
{
    unsigned counter = 0;
    for(unsigned i = 0; i < numPhases(); ++i, counter += phase(i).numSpecies())
        if(counter > index) return i;
    return numPhases();
}

auto ChemicalSystem::indicesElements(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, elements());
}

auto ChemicalSystem::indicesSpecies(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, species());
}

auto ChemicalSystem::indicesPhases(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, phases());
}

auto ChemicalSystem::indicesElementsInSpecies(Index index) const -> Indices
{
    Indices indices;
    for(Element element : species(index).elements())
        indices.push_back(indexElement(element.name()));
    return indices;
}

auto ChemicalSystem::indicesElementsInSpecies(const Indices& ispecies) const -> Indices
{
    std::set<Index> ielements;
    for(const Index& i : ispecies)
    {
        const Indices& indices = indicesElementsInSpecies(i);
        ielements.insert(indices.begin(), indices.end());
    }
    return Indices(ielements.begin(), ielements.end());
}

auto ChemicalSystem::indicesPhasesWithSpecies(const Indices& ispecies) const -> Indices
{
    Indices iphases;
    iphases.reserve(numPhases());
    for(Index i : ispecies)
        iphases.push_back(indexPhaseWithSpecies(i));
    return iphases;
}

auto ChemicalSystem::offset(Index iphase) const -> unsigned
{
    unsigned counter = 0;
    for(unsigned i = 0; i < iphase-1; ++i)
        counter += phase(i).species().size();
    return counter;
}

auto ChemicalSystem::gibbsEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.gibbs_energies(T, P);
}

auto ChemicalSystem::enthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->data.enthalpies(T, P);
}

auto ChemicalSystem::helmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.helmholtz_energies(T, P);
}

auto ChemicalSystem::entropies(double T, double P) const -> ThermoVector
{
    return pimpl->data.entropies(T, P);
}

auto ChemicalSystem::volumes(double T, double P) const -> ThermoVector
{
    return pimpl->data.volumes(T, P);
}

auto ChemicalSystem::internalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.internal_energies(T, P);
}

auto ChemicalSystem::heatCapacitiesCp(double T, double P) const -> ThermoVector
{
    return pimpl->data.heat_capacities_cp(T, P);
}

auto ChemicalSystem::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.concentrations(T, P, n);
}

auto ChemicalSystem::lnActivityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.ln_activity_coefficients(T, P, n);
}

auto ChemicalSystem::lnActivities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.ln_activities(T, P, n);
}

auto ChemicalSystem::chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.chemical_potentials(T, P, n);
}

auto ChemicalSystem::densities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.densities(T, P, n);
}

auto ChemicalSystem::b(const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    return W * n;
}

auto ChemicalSystem::bphase(const Vector& n, Index iphase) const -> Vector
{
    const Matrix& W = formulaMatrix();
    const unsigned first = offset(iphase);
    const unsigned size = phase(iphase).numSpecies();
    const auto Wp = cols(W, first, size);
    const auto np = cols(n, first, size);
    return Wp * np;
}

auto ChemicalSystem::bspecies(const Vector& n, const Indices& ispecies) const -> Vector
{
    const Matrix& W = formulaMatrix();
    Vector b = zeros(numElements());
    for(unsigned j = 0; j < b.rows(); ++j)
        for(Index i : ispecies)
            b[j] += W(j, i) * n[i];
    return b;
}

auto ChemicalSystem::nphase(const Vector& n, Index iphase) const -> Vector
{
    const unsigned first = offset(iphase);
    const unsigned size = phase(iphase).numSpecies();
    return cols(n, first, size);
}

} // namespace Reaktor
