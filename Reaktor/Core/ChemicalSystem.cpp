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

auto formulaMatrix(const std::vector<Element>& elements, const std::vector<Species>& species) -> Matrix
{
    const auto& num_elements = elements.size();
    const auto& num_species = species.size();
    Matrix W(num_elements, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned j = 0; j < num_elements; ++j)
            W(j, i) = species[i].elementAtoms(elements[j].name());
    return W;
}

auto raiseErrorIfThereAreSpeciesWithSameNames(const std::vector<Species>& species) -> void
{
    std::set<std::string> names;
    for(const Species& s : species)
    {
        if(names.count(s.name()))
            RuntimeError("Cannot initialize the ChemicalSystem instance.",
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
    /// The list of phases in the system
    std::vector<Phase> phases;

    /// The list of species in the system
    std::vector<Species> species;

    /// The list of elements in the system
    std::vector<Element> elements;

    /// The model configuration of the system
    ChemicalSystemModel model;

    /// The formula matrix of the system
    Matrix formula_matrix;

    Impl()
    {}

    Impl(const std::vector<Phase>& _phases)
    : phases(_phases), species(collectSpecies(phases)), elements(collectElements(species))
    {
        raiseErrorIfThereAreSpeciesWithSameNames(species);

        formula_matrix = Reaktor::formulaMatrix(elements, species);

        model.standard_gibbs_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardGibbsEnergy(T, P);
            return res;
        };

        model.standard_helmholtz_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHelmholtzEnergy(T, P);
            return res;
        };

        model.standard_internal_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardInternalEnergy(T, P);
            return res;
        };

        model.standard_enthalpy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEnthalpy(T, P);
            return res;
        };

        model.standard_entropy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEntropy(T, P);
            return res;
        };

        model.standard_volume_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardVolume(T, P);
            return res;
        };

        model.standard_heat_capacity_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHeatCapacity(T, P);
            return res;
        };

        model.concentration_fn = [&](double T, double P, const Vector& n)
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

        model.activity_coefficient_fn = [&](double T, double P, const Vector& n)
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

        model.activity_fn = [&](double T, double P, const Vector& n)
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

        model.chemical_potential_fn = [&](double T, double P, const Vector& n)
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

        model.phase_molar_volume_fn = [&](double T, double P, const Vector& n)
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
    }

    Impl(const std::vector<Phase>& _phases, const ChemicalSystemModel& _model)
    : Impl(_phases)
    {
        if(_model.standard_gibbs_energy_fn)
            model.standard_gibbs_energy_fn = _model.standard_gibbs_energy_fn;

        if(_model.standard_helmholtz_energy_fn)
            model.standard_helmholtz_energy_fn = _model.standard_helmholtz_energy_fn;

        if(_model.standard_internal_energy_fn)
            model.standard_internal_energy_fn = _model.standard_internal_energy_fn;

        if(_model.standard_enthalpy_fn)
            model.standard_enthalpy_fn = _model.standard_enthalpy_fn;

        if(_model.standard_entropy_fn)
            model.standard_entropy_fn = _model.standard_entropy_fn;

        if(_model.standard_volume_fn)
            model.standard_volume_fn = _model.standard_volume_fn;

        if(_model.standard_heat_capacity_fn)
            model.standard_heat_capacity_fn = _model.standard_heat_capacity_fn;

        if(_model.concentration_fn)
            model.concentration_fn = _model.concentration_fn;

        if(_model.activity_coefficient_fn)
            model.activity_coefficient_fn = _model.activity_coefficient_fn;

        if(_model.activity_fn)
            model.activity_fn = _model.activity_fn;

        if(_model.chemical_potential_fn)
            model.chemical_potential_fn = _model.chemical_potential_fn;

        if(_model.phase_molar_volume_fn)
            model.phase_molar_volume_fn = _model.phase_molar_volume_fn;
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases, const ChemicalSystemModel& model)
: pimpl(new Impl(phases, model))
{}

ChemicalSystem::~ChemicalSystem()
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

auto ChemicalSystem::indexElement(std::string name) const -> Index
{
    return index(name, elements());
}

auto ChemicalSystem::indexElementWithError(std::string name) const -> Index
{
    const Index index = indexElement(name);

    Assert(index < numElements(),
        "Cannot get the index of the element `" + name + "`.",
        "There is no element called `" + name + "` in the system.");

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
        "There is no element called `" + name + "` in the system.");

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
        "There is no phase called `" + name + "` in the system.");

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
    return pimpl->model.standard_gibbs_energy_fn(T, P);
}

auto ChemicalSystem::standardEnthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_helmholtz_energy_fn(T, P);
}

auto ChemicalSystem::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_internal_energy_fn(T, P);
}

auto ChemicalSystem::standardEntropies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_enthalpy_fn(T, P);
}

auto ChemicalSystem::standardVolumes(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_entropy_fn(T, P);
}

auto ChemicalSystem::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_volume_fn(T, P);
}

auto ChemicalSystem::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    return pimpl->model.standard_heat_capacity_fn(T, P);
}

auto ChemicalSystem::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->model.concentration_fn(T, P, n);
}

auto ChemicalSystem::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->model.activity_coefficient_fn(T, P, n);
}

auto ChemicalSystem::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->model.activity_fn(T, P, n);
}

auto ChemicalSystem::chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->model.chemical_potential_fn(T, P, n);
}

auto ChemicalSystem::phaseMolarVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->model.phase_molar_volume_fn(T, P, n);
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
