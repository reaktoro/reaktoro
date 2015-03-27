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

#include "Multiphase.hpp"

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

auto raiseErrorIfThereAreSpeciesWithSameNames(const std::vector<Species>& species) -> void
{
    std::set<std::string> names;
    for(const Species& s : species)
    {
        if(names.count(s.name()))
            RuntimeError("Cannot initialize the Multiphase instance.",
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

struct Multiphase::Impl
{
    /// The list of phases in the multiphase system
    std::vector<Phase> phases;

    /// The list of species in the multiphase system
    std::vector<Species> species;

    /// The list of elements in the multiphase system
    std::vector<Element> elements;

    /// The formula matrix of the multiphase system
    Matrix formula_matrix;

    /// The function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_gibbs_energy_fn;

    /// The function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_helmholtz_energy_fn;

    /// The function for the apparent standard molar internal energies of the species (in units of J/mol).
    ThermoVectorFunction standard_internal_energy_fn;

    /// The function for the apparent standard molar enthalpies of the species (in units of J/mol).
    ThermoVectorFunction standard_enthalpy_fn;

    /// The function for the standard molar entropies of the species (in units of J/K).
    ThermoVectorFunction standard_entropy_fn;

    /// The function for the standard molar volumes of the species (in units of m3/mol).
    ThermoVectorFunction standard_volume_fn;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoVectorFunction standard_heat_capacity_fn;

    /// The function for the concentrations of the species (no uniform units).
    ChemicalVectorFunction concentration_fn;

    /// The function for the natural log of the activity coefficients of the species.
    ChemicalVectorFunction activity_coefficient_fn;

    /// The function for the natural log of the activities of the species.
    ChemicalVectorFunction activity_fn;

    /// The function for the chemical potentials of the species (in units of J/mol).
    ChemicalVectorFunction chemical_potential_fn;

    /// The function for the molar volumes of the phases (in units of m3/mol).
    ChemicalVectorFunction phase_molar_volume_fn;

    Impl()
    {}

    Impl(const std::vector<Phase>& _phases)
    : phases(_phases), species(collectSpecies(phases)), elements(collectElements(species))
    {
        raiseErrorIfThereAreSpeciesWithSameNames(species);

        formula_matrix = Reaktor::formulaMatrix(elements, species);

        standard_gibbs_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardGibbsEnergy(T, P);
            return res;
        };

        standard_helmholtz_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHelmholtzEnergy(T, P);
            return res;
        };

        standard_internal_energy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardInternalEnergy(T, P);
            return res;
        };

        standard_enthalpy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEnthalpy(T, P);
            return res;
        };

        standard_entropy_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardEntropy(T, P);
            return res;
        };

        standard_volume_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardVolume(T, P);
            return res;
        };

        standard_heat_capacity_fn = [&](double T, double P)
        {
            ThermoVector res(species.size());
            for(unsigned i = 0; i < species.size(); ++i)
                res.row(i) = species[i].standardHeatCapacity(T, P);
            return res;
        };

        concentration_fn = [&](double T, double P, const Vector& n)
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

        activity_coefficient_fn = [&](double T, double P, const Vector& n)
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

        activity_fn = [&](double T, double P, const Vector& n)
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

        chemical_potential_fn = [&](double T, double P, const Vector& n)
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

        phase_molar_volume_fn = [&](double T, double P, const Vector& n)
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
};

Multiphase::Multiphase()
: pimpl(new Impl())
{}

Multiphase::Multiphase(const Multiphase& other)
: pimpl(new Impl(*other.pimpl))
{}

Multiphase::Multiphase(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

Multiphase::~Multiphase()
{}

auto Multiphase::operator=(Multiphase other) -> Multiphase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Multiphase::setStandardGibbsEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_gibbs_energy_fn = function;
}

auto Multiphase::setStandardEnthalpyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_helmholtz_energy_fn = function;
}

auto Multiphase::setStandardHelmholtzEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_internal_energy_fn = function;
}

auto Multiphase::setStandardInternalEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_enthalpy_fn = function;
}

auto Multiphase::setStandardEntropyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_entropy_fn = function;
}

auto Multiphase::setStandardVolumeFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_volume_fn = function;
}

auto Multiphase::setStandardHeatCapacityFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_heat_capacity_fn = function;
}

auto Multiphase::setConcentrationFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->concentration_fn = function;
}

auto Multiphase::setActivityCoefficientFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_coefficient_fn = function;
}

auto Multiphase::setActivityFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_fn = function;
}

auto Multiphase::setChemicalPotentialFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->chemical_potential_fn = function;
}

auto Multiphase::setPhaseMolarVolumeFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->phase_molar_volume_fn = function;
}

auto Multiphase::numElements() const -> unsigned
{
    return elements().size();
}

auto Multiphase::numSpecies() const -> unsigned
{
    return species().size();
}

auto Multiphase::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return phase(iphase).numSpecies();
}

auto Multiphase::numPhases() const -> unsigned
{
    return phases().size();
}

auto Multiphase::element(Index index) const -> const Element&
{
    return elements()[index];
}

auto Multiphase::element(std::string name) const -> const Element&
{
    return element(indexElementWithError(name));
}

auto Multiphase::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto Multiphase::species(Index index) const -> const Species&
{
    return species()[index];
}

auto Multiphase::species(std::string name) const -> const Species&
{
    return species(indexSpeciesWithError(name));
}

auto Multiphase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Multiphase::phase(Index index) const -> const Phase&
{
    return phases()[index];
}

auto Multiphase::phase(std::string name) const -> const Phase&
{
    return phase(indexPhaseWithError(name));
}

auto Multiphase::phases() const -> const std::vector<Phase>&
{
    return pimpl->phases;
}

auto Multiphase::formulaMatrix() const -> const Matrix&
{
    return pimpl->formula_matrix;
}

auto Multiphase::indexElement(std::string name) const -> Index
{
    return index(name, elements());
}

auto Multiphase::indexElementWithError(std::string name) const -> Index
{
    const Index index = indexElement(name);

    Assert(index < numElements(),
        "Cannot get the index of the element `" + name + "`.",
        "There is no element called `" + name + "` in the system.");

    return index;
}

auto Multiphase::indexSpecies(std::string name) const -> Index
{
    return index(name, species());
}

auto Multiphase::indexSpeciesWithError(std::string name) const -> Index
{
    const Index index = indexSpecies(name);

    Assert(index < numSpecies(),
        "Cannot get the index of the species `" + name + "`.",
        "There is no element called `" + name + "` in the system.");

    return index;
}

auto Multiphase::indexPhase(std::string name) const -> Index
{
    return index(name, phases());
}

auto Multiphase::indexPhaseWithError(std::string name) const -> Index
{
    const Index index = indexPhase(name);

    Assert(index < numPhases(),
        "Cannot get the index of the phase `" + name + "`.",
        "There is no phase called `" + name + "` in the system.");

    return index;
}

auto Multiphase::indexPhaseWithSpecies(Index index) const -> Index
{
    unsigned counter = 0;
    for(unsigned i = 0; i < numPhases(); ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > index) return i;
    }
    return numPhases();
}

auto Multiphase::indicesElements(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, elements());
}

auto Multiphase::indicesSpecies(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, species());
}

auto Multiphase::indicesPhases(const std::vector<std::string>& names) const -> Indices
{
    return indices(names, phases());
}

auto Multiphase::indicesElementsInSpecies(Index index) const -> Indices
{
    Indices indices;
    for(const auto& pair : species(index).elements())
        indices.push_back(indexElement(pair.first.name()));
    return indices;
}

auto Multiphase::indicesElementsInSpecies(const Indices& ispecies) const -> Indices
{
    std::set<Index> ielements;
    for(const Index& i : ispecies)
    {
        const Indices& indices = indicesElementsInSpecies(i);
        ielements.insert(indices.begin(), indices.end());
    }
    return Indices(ielements.begin(), ielements.end());
}

auto Multiphase::indicesPhasesWithSpecies(const Indices& ispecies) const -> Indices
{
    Indices iphases;
    iphases.reserve(numPhases());
    for(Index i : ispecies)
        iphases.push_back(indexPhaseWithSpecies(i));
    return iphases;
}

auto Multiphase::indexFirstSpeciesInPhase(Index iphase) const -> unsigned
{
    unsigned counter = 0;
    for(unsigned i = 0; i < iphase; ++i)
        counter += phase(i).species().size();
    return counter;
}

auto Multiphase::standardGibbsEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->standard_gibbs_energy_fn(T, P);
}

auto Multiphase::standardEnthalpies(double T, double P) const -> ThermoVector
{
    return pimpl->standard_helmholtz_energy_fn(T, P);
}

auto Multiphase::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->standard_internal_energy_fn(T, P);
}

auto Multiphase::standardEntropies(double T, double P) const -> ThermoVector
{
    return pimpl->standard_enthalpy_fn(T, P);
}

auto Multiphase::standardVolumes(double T, double P) const -> ThermoVector
{
    return pimpl->standard_entropy_fn(T, P);
}

auto Multiphase::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    return pimpl->standard_volume_fn(T, P);
}

auto Multiphase::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    return pimpl->standard_heat_capacity_fn(T, P);
}

auto Multiphase::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->concentration_fn(T, P, n);
}

auto Multiphase::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->activity_coefficient_fn(T, P, n);
}

auto Multiphase::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->activity_fn(T, P, n);
}

auto Multiphase::chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->chemical_potential_fn(T, P, n);
}

auto Multiphase::phaseMolarVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->phase_molar_volume_fn(T, P, n);
}

auto Multiphase::phaseDensities(double T, double P, const Vector& n) const -> ChemicalVector
{
    // TODO Implement this
    RuntimeError("Multiphase::phaseDensities has not been implemented yet.", "");
}

auto Multiphase::phaseMolarAmounts(const Vector& n) const -> ChemicalVector
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

auto Multiphase::phaseMassAmounts(const Vector& n) const -> Vector
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

auto Multiphase::phaseVolumes(double T, double P, const Vector& n) const -> ChemicalVector
{
    const ChemicalVector vphases = phaseMolarVolumes(T, P, n);
    const ChemicalVector nphases = phaseMolarAmounts(n);
    return vphases % nphases;
}

auto Multiphase::elementAmounts(const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    return W * n;
}

auto Multiphase::elementAmountsInPhase(Index iphase, const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    const unsigned first = indexFirstSpeciesInPhase(iphase);
    const unsigned size = numSpeciesInPhase(iphase);
    const auto Wp = cols(W, first, size);
    const auto np = cols(n, first, size);
    return Wp * np;
}

auto Multiphase::elementAmountsInSpecies(const Indices& ispecies, const Vector& n) const -> Vector
{
    const Matrix& W = formulaMatrix();
    Vector b = zeros(W.rows());
    for(Index i : ispecies)
        b += W.col(i) * n[i];
    return b;
}

auto Multiphase::elementAmount(Index ielement, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    return W.row(ielement) * n;
}

auto Multiphase::elementAmountInPhase(Index ielement, Index iphase, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    const unsigned first = indexFirstSpeciesInPhase(iphase);
    const unsigned size = numSpeciesInPhase(iphase);
    const auto Wp = cols(W, first, size);
    const auto np = cols(n, first, size);
    return dot(Wp.row(ielement), np);
}

auto Multiphase::elementAmountInSpecies(Index ielement, const Indices& ispecies, const Vector& n) const -> double
{
    const Matrix& W = formulaMatrix();
    double bval = 0.0;
    for(Index i : ispecies)
        bval += W(ielement, i) * n[i];
    return bval;
}

auto operator<<(std::ostream& out, const Multiphase& multiphase) -> std::ostream&
{
    const auto& phases = multiphase.phases();
    for(unsigned i = 0; i < phases.size(); ++i)
    {
        out << "Phase(" << i << "): " << phases[i].name() << std::endl;
        const auto& species = phases[i].species();
        for(unsigned i = 0; i < species.size(); ++i)
        {
            const auto name = species[i].name();
            const auto idx  = multiphase.indexSpecies(name);
            out << std::setw(5) << std::left << idx;
            out << std::setw(30) << std::left << name;
            out << std::endl;
        }
    }

    out << std::endl;

    out << "Elements:" << std::endl;
    const auto& elements = multiphase.elements();
    for(unsigned i = 0; i < elements.size(); ++i)
        out << (i > 0 ? ", " : "") << i << ":" << elements[i].name();

    out << std::endl;

    return out;
}

} // namespace Reaktor
