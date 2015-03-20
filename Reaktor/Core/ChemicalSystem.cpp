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
#include <iostream>
#include <iomanip>
#include <set>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/OptimizationUtils.hpp>
#include <Reaktor/Core/CoreUtils.hpp>

namespace Reaktor {
namespace {

auto formulaMatrix(
    const std::vector<Element>& elements,
    const std::vector<Species>& species) -> Matrix
{
    const auto& num_elements = elements.size();
    const auto& num_species = species.size();
    Matrix W(num_elements, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned j = 0; j < num_elements; ++j)
            W(j, i) = atoms(elements[j], species[i]);
    return W;
}

auto connectivity(
    const std::vector<Element>& elements,
    const std::vector<Species>& species,
    const std::vector<Phase>& phases) -> Connectivity
{
    const unsigned num_elements = elements.size();
    const unsigned num_species = species.size();
    const unsigned num_phases = phases.size();

    Connectivity c;

    c.element_to_species.resize(num_elements);
    c.species_to_elements.resize(num_species);
    for(unsigned j = 0; j < num_elements; ++j)
        for(unsigned i = 0; i < num_species; ++i)
            if(species[i].elements().count(elements[j])) {
                c.element_to_species[j].push_back(i);
                c.species_to_elements[i].push_back(j); }

    c.species_to_phase.resize(num_species);
    c.phase_to_species.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned i = 0; i < num_species; ++i)
            if(contains(species[i], phases[k].species())) {
                c.species_to_phase[i] = k;
                c.phase_to_species[k].push_back(i);
            }

    c.element_to_phases.resize(num_elements);
    c.phase_to_elements.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned j = 0; j < num_elements; ++j)
            if(contains(elements[j], phases[k].elements())) {
                c.element_to_phases[j].push_back(k);
                c.phase_to_elements[k].push_back(j);
            }

    return c;
}

auto optimize(const ChemicalSystemData& data) -> ChemicalSystemData
{
    // TODO Optimize other functions as well
    ChemicalSystemData optimized = data;
    optimized.ln_activity_coefficients = memoizeLast(data.ln_activity_coefficients);
    optimized.ln_activities = memoizeLast(data.ln_activities);
    optimized.chemical_potentials = memoizeLast(data.chemical_potentials);
    optimized.phase_molar_volumes = memoizeLast(data.phase_molar_volumes);
    return optimized;
}

auto checkSpeciesWithSameNames(const std::vector<Species>& species) -> void
{
    std::set<std::string> names;
    for(const Species& s : species)
    {
        if(names.count(s.name()))
            RuntimeError("Cannot initialise the ChemicalSystem instance.",
                "The species `" + s.name() + "` has more than one occurrence.");
        names.insert(s.name());
    }
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

    /// The connectivity of the chemical system
    Connectivity connectivity;

    Impl()
    {}

    Impl(const ChemicalSystemData& data)
    : data(optimize(data)), species(collectSpecies(data.phases)), elements(collectElements(species))
    {
        checkSpeciesWithSameNames(species);

        formula_matrix = Reaktor::formulaMatrix(elements, species);
        connectivity = Reaktor::connectivity(elements, species, data.phases);
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

auto ChemicalSystem::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return phase(iphase).numSpecies();
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
    return element(indexElementWithError(name));
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
    return species(indexSpeciesWithError(name));
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
    return phase(indexPhaseWithError(name));
}

auto ChemicalSystem::phases() const -> const std::vector<Phase>&
{
    return pimpl->data.phases;
}

auto ChemicalSystem::formulaMatrix() const -> const Matrix&
{
    return pimpl->formula_matrix;
}

auto ChemicalSystem::connectivity() const -> const Connectivity&
{
    return pimpl->connectivity;
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
    for(unsigned i = 0; i < numPhases(); ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > index) return i;
    }
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
    for(const auto& pair : species(index).elements())
        indices.push_back(indexElement(pair.first.name()));
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

auto ChemicalSystem::indexFirstSpeciesInPhase(Index iphase) const -> unsigned
{
    unsigned counter = 0;
    for(unsigned i = 0; i < iphase; ++i)
        counter += phase(i).species().size();
    return counter;
}

auto ChemicalSystem::standardGibbsEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_gibbs_energies(T, P);
}

auto ChemicalSystem::standardEnthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_enthalpies(T, P);
}

auto ChemicalSystem::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_helmholtz_energies(T, P);
}

auto ChemicalSystem::standardEntropies(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_entropies(T, P);
}

auto ChemicalSystem::standardVolumes(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_volumes(T, P);
}

auto ChemicalSystem::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_internal_energies(T, P);
}

auto ChemicalSystem::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    return pimpl->data.standard_heat_capacities(T, P);
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

auto ChemicalSystem::phaseMolarVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->data.phase_molar_volumes(T, P, n);
}

auto ChemicalSystem::phaseDensities(double T, double P, const Vector& n) const -> ChemicalVector
{
    RuntimeError("Cannot calculate phase densities.", "Method ChemicalSystem::phaseDensities has not been implemented yet.");
    return ChemicalVector();
}

auto ChemicalSystem::phaseTotalAmounts(const Vector& n) const -> Vector
{
    const unsigned num_phases = numPhases();
    Vector nphases(num_phases);
    unsigned offset = 0;
    for(unsigned i = 0; i < num_phases; ++i)
    {
        const unsigned size = numSpeciesInPhase(i);
        const auto np = rows(n, offset, size);
        nphases[i] = sum(np);
        offset += size;
    }
    return nphases;
}

auto ChemicalSystem::phaseVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    const unsigned num_species = numSpecies();
    const unsigned num_phases = numPhases();
    const Vector vphases = phaseMolarVolumes(T, P, n).val;
    const Vector nphases = phaseTotalAmounts(n);
    const Vector Vphases = vphases % nphases;
    const Vector zero_vec = zeros(num_phases);
    const Matrix zero_mat = zeros(num_phases, num_species);
    return ChemicalVector(Vphases, zero_vec, zero_vec, zero_mat);
}

auto ChemicalSystem::elementAmounts(const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    return W * n;
}

auto ChemicalSystem::elementAmountsInPhase(Index iphase, const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    const unsigned first = indexFirstSpeciesInPhase(iphase);
    const unsigned size = numSpeciesInPhase(iphase);
    const auto Wp = cols(W, first, size);
    const auto np = cols(n, first, size);
    return Wp * np;
}

auto ChemicalSystem::elementAmountsInSpecies(const Indices& ispecies, const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    Vector b = zeros(W.rows());
    for(Index i : ispecies)
        b += W.col(i) * n[i];
    return b;
}

auto ChemicalSystem::elementAmount(Index ielement, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    return W.row(ielement) * n;
}

auto ChemicalSystem::elementAmountInPhase(Index ielement, Index iphase, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    const unsigned first = indexFirstSpeciesInPhase(iphase);
    const unsigned size = numSpeciesInPhase(iphase);
    const auto Wp = cols(W, first, size);
    const auto np = cols(n, first, size);
    return dot(Wp.row(ielement), np);
}

auto ChemicalSystem::elementAmountInSpecies(Index ielement, const Indices& ispecies, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    double bval = 0.0;
    for(Index i : ispecies)
        bval += W(ielement, i) * n[i];
    return bval;
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    const auto& phases = system.phases();
    for(unsigned i = 0; i < phases.size(); ++i)
    {
        out << "Phase(" << i << "): " << phases[i].name() << std::endl;
        const auto& species = phases[i].species();
        for(unsigned i = 0; i < species.size(); ++i)
        {
            const auto name = species[i].name();
            const auto idx  = system.indexSpecies(name);
            out << std::setw(5) << std::left << idx;
            out << std::setw(30) << std::left << name;
            out << std::endl;
        }
    }

    out << std::endl;

    out << "Elements:" << std::endl;
    const auto& elements = system.elements();
    for(unsigned i = 0; i < elements.size(); ++i)
        out << (i > 0 ? ", " : "") << i << ":" << elements[i].name();

    out << std::endl;

    return out;
}

} // namespace Reaktor
