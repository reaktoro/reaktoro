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
#include <Reaktor/Core/Utils.hpp>

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
            W(j, i) = species[i].elementAtoms(elements[j].name());
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

auto raiseErrorIfThereAreSpeciesWithSameNames(const std::vector<Species>& species) -> void
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

auto collectElements(const std::vector<Species>& species) -> std::vector<Element>
{
    std::set<Element> elements;
    for(const Species& iter : species)
        for(const auto& pair : iter.elements())
            elements.insert(pair.first);
    return std::vector<Element>(elements.begin(), elements.end());
}

auto collectSpecies(const std::vector<Phase>& phases) -> std::vector<Species>
{
    unsigned num_species = 0;
    for(const Phase& phase : phases)
        num_species += phase.species().size();

    std::vector<Species> list;
    list.reserve(num_species);
    for(const Phase& phase : phases)
        for(const Species& iter : phase.species())
            list.push_back(iter);
    return list;
}

} // namespace

struct ChemicalSystem::Impl
{
    /// The list of phases in the chemical system
    std::vector<Phase> phases;

    /// The list of species in the chemical system
    std::vector<Species> species;

    /// The list of elements in the chemical system
    std::vector<Element> elements;

    /// The thermodynamic and physical property functions of the chemical system
    ChemicalModels models;

    /// The formula matrix of the chemical system
    Matrix formula_matrix;

    /// The connectivity of the chemical system
    Connectivity connectivity;

    Impl()
    {}

    Impl(const std::vector<Phase>& _phases)
    : Impl(_phases, {})
    {
        ThermoVectorFunction standard_gibbs_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardGibbsEnergy(T, P);
            return res;
        };

        ThermoVectorFunction standard_helmholtz_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHelmholtzEnergy(T, P);
            return res;
        };

        ThermoVectorFunction standard_internal_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardInternalEnergy(T, P);
            return res;
        };

        ThermoVectorFunction standard_enthalpy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEnthalpy(T, P);
            return res;
        };

        ThermoVectorFunction standard_entropy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEntropy(T, P);
            return res;
        };

        ThermoVectorFunction standard_volume_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardVolume(T, P);
            return res;
        };

        ThermoVectorFunction standard_heat_capacity_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHeatCapacity(T, P);
            return res;
        };

        ChemicalVectorFunction concentration_fn = [&](double T, double P, const Vector& n)
        {
            ChemicalVector res(species.size(), species.size());
            unsigned offset = 0;
            for(unsigned i = 0; i < phases.size(); ++i)
            {
                const unsigned size = phases[i].numSpecies();
                const auto np = rows(n, offset, size);
                res.rows(offset, offset, size, size) = phases[i].concentrations(T, P, np);
                offset += size;
            }
            return res;
        };

        ChemicalVectorFunction activity_coefficient_fn = [&](double T, double P, const Vector& n)
        {
            ChemicalVector res(species.size(), species.size());
            unsigned offset = 0;
            for(unsigned i = 0; i < phases.size(); ++i)
            {
                const unsigned size = phases[i].numSpecies();
                const auto np = rows(n, offset, size);
                res.rows(offset, offset, size, size) = phases[i].activityCoefficients(T, P, np);
                offset += size;
            }
            return res;
        };

        ChemicalVectorFunction activity_fn = [&](double T, double P, const Vector& n)
        {
            ChemicalVector res(species.size(), species.size());
            unsigned offset = 0;
            for(unsigned i = 0; i < phases.size(); ++i)
            {
                const unsigned size = phases[i].numSpecies();
                const auto np = rows(n, offset, size);
                res.rows(offset, offset, size, size) = phases[i].activities(T, P, np);
                offset += size;
            }
            return res;
        };

        ChemicalVectorFunction chemical_potential_fn = [&](double T, double P, const Vector& n)
        {
            ChemicalVector res(species.size(), species.size());
            unsigned offset = 0;
            for(unsigned i = 0; i < phases.size(); ++i)
            {
                const unsigned size = phases[i].numSpecies();
                const auto np = rows(n, offset, size);
                res.rows(offset, offset, size, size) = phases[i].chemicalPotentials(T, P, np);
                offset += size;
            }
            return res;
        };

        ChemicalVectorFunction phase_molar_volume_fn = [&](double T, double P, const Vector& n)
        {
            ChemicalVector res(phases.size(), species.size());
            unsigned offset = 0;
            for(unsigned i = 0; i < phases.size(); ++i)
            {
                const unsigned size = phases[i].numSpecies();
                const auto np = rows(n, offset, size);
                res.row(i, offset, size) = phases[i].molarVolume(T, P, np);
                offset += size;
            }
            return res;
        };

        models.setStandardGibbsEnergyFunction(standard_gibbs_energy_fn);
        models.setStandardHelmholtzEnergyFunction(standard_helmholtz_energy_fn);
        models.setStandardInternalEnergyFunction(standard_internal_energy_fn);
        models.setStandardEnthalpyFunction(standard_enthalpy_fn);
        models.setStandardEntropyFunction(standard_entropy_fn);
        models.setStandardVolumeFunction(standard_volume_fn);
        models.setStandardHeatCapacityFunction(standard_heat_capacity_fn);
        models.setConcentrationFunction(concentration_fn);
        models.setActivityCoefficientFunction(activity_coefficient_fn);
        models.setActivityFunction(activity_fn);
        models.setChemicalPotentialFunction(chemical_potential_fn);
        models.setPhaseMolarVolumeFunction(phase_molar_volume_fn);
    }

    Impl(const std::vector<Phase>& _phases, const ChemicalModels& _models)
    : phases(_phases), species(collectSpecies(phases)),
      elements(collectElements(species)), models(_models)
    {
        raiseErrorIfThereAreSpeciesWithSameNames(species);
        formula_matrix = Reaktor::formulaMatrix(elements, species);
        connectivity = Reaktor::connectivity(elements, species, phases);
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases, const ChemicalModels& models)
: pimpl(new Impl(phases, models))
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
    return pimpl->phases;
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
        "Cannot get the index of the element `" + name + "`.",
        "There is no element called `" + name + "` in the chemical system.");

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
        "Cannot get the index of the species `" + name + "`.",
        "There is no element called `" + name + "` in the chemical system.");

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
        "Cannot get the index of the phase `" + name + "`.",
        "There is no phase called `" + name + "` in the chemical system.");

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
    return pimpl->models.standardGibbsEnergyFunction()(T, P);
}

auto ChemicalSystem::standardEnthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardEnthalpyFunction()(T, P);
}

auto ChemicalSystem::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardHelmholtzEnergyFunction()(T, P);
}

auto ChemicalSystem::standardEntropies(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardEntropyFunction()(T, P);
}

auto ChemicalSystem::standardVolumes(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardVolumeFunction()(T, P);
}

auto ChemicalSystem::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardInternalEnergyFunction()(T, P);
}

auto ChemicalSystem::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    return pimpl->models.standardHeatCapacityFunction()(T, P);
}

auto ChemicalSystem::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->models.concentrationFunction()(T, P, n);
}

auto ChemicalSystem::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->models.activityCoefficientFunction()(T, P, n);
}

auto ChemicalSystem::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->models.activityFunction()(T, P, n);
}

auto ChemicalSystem::chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->models.chemicalPotentialFunction()(T, P, n);
}

auto ChemicalSystem::phaseMolarVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->models.phaseMolarVolumeFunction()(T, P, n);
}

auto ChemicalSystem::phaseDensities(double T, double P, const Vector& n) const -> ChemicalVector
{
    // TODO Implement this
    RuntimeError("ChemicalSystem::phaseDensities has not been implemented yet.", "");
}

auto ChemicalSystem::phaseMolarAmounts(const Vector& n) const -> ChemicalVector
{
    const unsigned num_species = numSpecies();
    const unsigned num_phases = numPhases();
    ChemicalVector nphases(num_phases, num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < num_phases; ++i)
    {
        const unsigned size = numSpeciesInPhase(i);
        const auto np = rows(n, offset, size);
        nphases.val[i] = sum(np);
        nphases.ddn.row(i).segment(offset, size).fill(1.0);
        offset += size;
    }
    return nphases;
}

auto ChemicalSystem::phaseMassAmounts(const Vector& n) const -> Vector
{
    const unsigned num_phases = numPhases();
    const Vector m = molarMasses(species());
    Vector mass_phases(num_phases);
    unsigned offset = 0;
    for(unsigned i = 0; i < num_phases; ++i)
    {
        const unsigned size = numSpeciesInPhase(i);
        const auto np = rows(n, offset, size);
        const auto mp = rows(m, offset, size);
        mass_phases[i] = sum(np % mp);
        offset += size;
    }
    return mass_phases;
}

auto ChemicalSystem::phaseVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    const ChemicalVector vphases = phaseMolarVolumes(T, P, n);
    const ChemicalVector nphases = phaseMolarAmounts(n);
    return vphases % nphases;
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
