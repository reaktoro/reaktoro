// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

auto formulaMatrix(const std::vector<Element>& elements, const std::vector<Species>& species) -> Matrix
{
    const auto& num_elements = elements.size();
    const auto& num_species = species.size();
    Matrix W = zeros(num_elements, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned j = 0; j < num_elements; ++j)
            W(j, i) = species[i].elementCoefficient(elements[j].name());
    return W;
}

auto raiseErrorIfThereAreSpeciesWithSameNames(const std::vector<Species>& species) -> void
{
    std::set<std::string> names;
    for(const Species& s : species)
    {
        if(names.count(s.name()))
            RuntimeError("Could not initialize the ChemicalSystem instance.",
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

auto createThermoModel(const std::vector<Phase>& phases) -> ThermoModel
{
    // The number of phases
    unsigned nphases = phases.size();

    // The number of species in all phases
    unsigned nspecies = 0;
    for(unsigned i = 0; i < nphases; ++i)
        nspecies += phases[i].numSpecies();

    // Define the thermodynamic model function for the system
    ThermoModel model = [=](double T, double P) -> ThermoModelResult
    {
        // The result of this function initialized for `nspecies`
        ThermoModelResult res(nphases);

        // The offset index of the first species in each phase
        unsigned offset = 0;

        // Iterate over all phases and calculate their thermodynamic properties
        for(unsigned i = 0; i < nphases; ++i)
        {
            // The number of species in the current phase
            const unsigned size = phases[i].numSpecies();

            // Calculate the standard thermodynamic properties of the current phase
            PhaseThermoModelResult phase_res = phases[i].thermoModel()(T, P);

            // Set the standard thermodynamic properties of the species in the current phase
            res[i].standard_partial_molar_gibbs_energies     = phase_res.standard_partial_molar_gibbs_energies;
            res[i].standard_partial_molar_enthalpies         = phase_res.standard_partial_molar_enthalpies;
            res[i].standard_partial_molar_volumes            = phase_res.standard_partial_molar_volumes;
            res[i].standard_partial_molar_heat_capacities_cp = phase_res.standard_partial_molar_heat_capacities_cp;
            res[i].standard_partial_molar_heat_capacities_cv = phase_res.standard_partial_molar_heat_capacities_cv;

            // Update the index of the first species in the next phase
            offset += size;
        }

        return res;
    };

    return model;
}

auto createChemicalModel(const std::vector<Phase>& phases) -> ChemicalModel
{
    // The number of phases
    unsigned nphases = phases.size();

    // The number of species in all phases
    unsigned nspecies = 0;
    for(unsigned i = 0; i < nphases; ++i)
        nspecies += phases[i].numSpecies();

    // Define the chemical model function for the system
    ChemicalModel model = [=](double T, double P, const Vector& n) -> ChemicalModelResult
    {
        // The result of this function initialized for `nspecies` and `nphases`
        ChemicalModelResult res(nphases);

        // The offset index of the first species in each phase
        unsigned offset = 0;

        // Iterate over all phases and calculate their chemical properties
        for(unsigned i = 0; i < nphases; ++i)
        {
            // The number of species in the current phase
            const unsigned size = phases[i].numSpecies();

            // The molar amounts of the species in the current phase
            const Vector np = rows(n, offset, size);

            // Calculate the chemical properties of the current phase
            PhaseChemicalModelResult phase_res = phases[i].chemicalModel()(T, P, np);

            // Set the chemical properties of the current phase
            res[i].ln_activity_coefficients        = phase_res.ln_activity_coefficients;
            res[i].ln_activities                   = phase_res.ln_activities;
            res[i].molar_volume                    = phase_res.molar_volume;
            res[i].residual_molar_gibbs_energy     = phase_res.residual_molar_gibbs_energy;
            res[i].residual_molar_enthalpy         = phase_res.residual_molar_enthalpy;
            res[i].residual_molar_heat_capacity_cp = phase_res.residual_molar_heat_capacity_cp;
            res[i].residual_molar_heat_capacity_cv = phase_res.residual_molar_heat_capacity_cv;

            // Update the index of the first species in the next phase
            offset += size;
        }

        return res;
    };

    return model;
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

    /// The formula matrix of the system
    Matrix formula_matrix;

    // The molar masses of the species in each phase
    std::vector<Vector> molar_masses;

    /// The function that calculates the thermodynamic properties of the species in the system
    ThermoModel thermo_model;

    /// The function that calculates the chemical properties of the species and phases in the system
    ChemicalModel chemical_model;

    Impl()
    {}

    Impl(const std::vector<Phase>& phaselist)
    : phases(phaselist), species(collectSpecies(phases)), elements(collectElements(species))
    {
        // Check if there are species with same names
        raiseErrorIfThereAreSpeciesWithSameNames(species);

        // Initialize the default thermodynamic model function of the system
        thermo_model = createThermoModel(phases);

        // Initialize the default chemical model function of the system
        chemical_model = createChemicalModel(phases);

        // Initialize the formula matrix of the system
        formula_matrix = Reaktoro::formulaMatrix(elements, species);

        // Initialize the vector of molar masses of the species in each phase
        molar_masses.reserve(phases.size());
        for(const Phase& phase : phases)
            molar_masses.push_back(molarMasses(phase.species()));
    }

    /// Calculate the standard thermodynamic properties of the species
    auto properties(double T, double P) const -> ThermoProperties
    {
        // The number of species and phases in the system
        const unsigned nspecies = species.size();
        const unsigned nphases = phases.size();

        // Evaluate the custom thermodynamic model function
        ThermoModelResult res = thermo_model(T, P);

        // The standard thermodynamic properties of the species at (T, P)
        ThermoProperties prop(nspecies);

        // Set temperature and pressure
        prop.T = ThermoScalar::Temperature(T);
        prop.P = ThermoScalar::Pressure(P);

        // The offset index of the first species in each phase
        unsigned offset = 0;

        // Iterate over all phases and calculate their chemical properties
        for(unsigned i = 0; i < nphases; ++i)
        {
            // The number of species in the current phase
            const unsigned size = phases[i].numSpecies();

            // Set the standard thermodynamic properties of the species in the current phase
            prop.standard_partial_molar_gibbs_energies.rows(offset, size)     = res[i].standard_partial_molar_gibbs_energies;
            prop.standard_partial_molar_enthalpies.rows(offset, size)         = res[i].standard_partial_molar_enthalpies;
            prop.standard_partial_molar_volumes.rows(offset, size)            = res[i].standard_partial_molar_volumes;
            prop.standard_partial_molar_heat_capacities_cp.rows(offset, size) = res[i].standard_partial_molar_heat_capacities_cp;
            prop.standard_partial_molar_heat_capacities_cv.rows(offset, size) = res[i].standard_partial_molar_heat_capacities_cv;

            // Update the index of the first species in the next phase
            offset += size;
        }

        return prop;
    }

    /// Calculate the chemical and thermodynamic properties of the chemical system
    auto properties(double T, double P, const Vector& n) const -> ChemicalProperties
    {
        // The number of phases and species in the system
        const unsigned nphases = phases.size();
        const unsigned nspecies = species.size();

        // The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
        ChemicalProperties prop(nspecies, nphases);

        // Set temperature, pressure and composition
        prop.T = ThermoScalar::Temperature(T);
        prop.P = ThermoScalar::Pressure(P);
        prop.n = ChemicalVector::Composition(n);

        // Calculate the thermodynamic properties of the system
        ThermoModelResult tres = thermo_model(T, P);

        // Calculate the chemical properties of the system
        ChemicalModelResult cres = chemical_model(T, P, n);

        // The offset index of the first species in each phase
        unsigned offset = 0;

        // Iterate over all phases and calculate their chemical properties
        for(unsigned i = 0; i < nphases; ++i)
        {
            // The number of species in the current phase
            const unsigned size = phases[i].numSpecies();

            // The molar composition of the species in the current phase
            const ChemicalVector np = ChemicalVector::Composition(rows(n, offset, size));

            // The molar fractions of the species in the current phase
            const ChemicalVector xp = molarFractions(np);

            // Set the standard thermodynamic properties of the species in the current phase
            prop.standard_partial_molar_gibbs_energies.rows(offset, size)     = tres[i].standard_partial_molar_gibbs_energies;
            prop.standard_partial_molar_enthalpies.rows(offset, size)         = tres[i].standard_partial_molar_enthalpies;
            prop.standard_partial_molar_volumes.rows(offset, size)            = tres[i].standard_partial_molar_volumes;
            prop.standard_partial_molar_heat_capacities_cp.rows(offset, size) = tres[i].standard_partial_molar_heat_capacities_cp;
            prop.standard_partial_molar_heat_capacities_cv.rows(offset, size) = tres[i].standard_partial_molar_heat_capacities_cv;

            // Set the molar fractions of the species in the current phase
            prop.molar_fractions.rows(offset, offset, size, size) = xp;

            // Set the ln activities and ln activity coefficients of the species in the current phase
            prop.ln_activity_coefficients.rows(offset, offset, size, size) = cres[i].ln_activity_coefficients;
            prop.ln_activities.rows(offset, offset, size, size) = cres[i].ln_activities;

            // Calculate the ideal contribution for the thermodynamic properties of the phase
            prop.phase_molar_gibbs_energies.row(i, offset, size)     = sum(xp % tres[i].standard_partial_molar_gibbs_energies);
            prop.phase_molar_enthalpies.row(i, offset, size)         = sum(xp % tres[i].standard_partial_molar_enthalpies);
            prop.phase_molar_volumes.row(i, offset, size)            = sum(xp % tres[i].standard_partial_molar_volumes);
            prop.phase_molar_heat_capacities_cp.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_heat_capacities_cp);
            prop.phase_molar_heat_capacities_cv.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_heat_capacities_cv);

            // Check if the molar volume of the phase was calculated
            if(cres[i].molar_volume.val > 0.0)
                prop.phase_molar_volumes.row(i, offset, size) = cres[i].molar_volume;

            // Add the non-ideal residual contribution to the thermodynamic properties of the phase
            prop.phase_molar_gibbs_energies.row(i, offset, size)     += cres[i].residual_molar_gibbs_energy;
            prop.phase_molar_enthalpies.row(i, offset, size)         += cres[i].residual_molar_enthalpy;
            prop.phase_molar_heat_capacities_cp.row(i, offset, size) += cres[i].residual_molar_heat_capacity_cp;
            prop.phase_molar_heat_capacities_cv.row(i, offset, size) += cres[i].residual_molar_heat_capacity_cv;

            // Set the molar amount and mass of the current phase
            prop.phase_moles.row(i, offset, size) = sum(np);
            prop.phase_masses.row(i, offset, size) = sum(molar_masses[i] % np);

            // Update the index of the first species in the next phase
            offset += size;
        }

        return prop;
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

ChemicalSystem::~ChemicalSystem()
{}

auto ChemicalSystem::setThermoModel(const ThermoModel& model) -> void
{
    pimpl->thermo_model = model;
}

auto ChemicalSystem::setChemicalModel(const ChemicalModel& model) -> void
{
    pimpl->chemical_model = model;
}

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
    Assert(index < numPhases(),
        "Could not get a reference to a Phase instance with given index.",
        "The given index " + std::to_string(index) + " is out of bounds.")
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
        "Could not get the index of element `" + name + "`.",
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
        "Could not get the index of species `" + name + "`.",
        "There is no species called `" + name + "` in the system.");
    return index;
}

auto ChemicalSystem::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, species());
}

auto ChemicalSystem::indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index
{
    const Index index = indexSpeciesAny(names);
    Assert(index < numSpecies(),
        "Could not get the index of the species with "
        "any of the following names `" + join(names, ", ") + "`.",
        "There is no species in the system with any of these names.");
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
        "Could not get the index of phase `" + name + "`.",
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
    const auto np = rows(n, first, size);
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
    const auto np = rows(n, first, size);
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

auto ChemicalSystem::properties(double T, double P) const -> ThermoProperties
{
    return pimpl->properties(T, P);
}

auto ChemicalSystem::properties(double T, double P, const Vector& n) const -> ChemicalProperties
{
    return pimpl->properties(T, P, n);
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    const auto& phases = system.phases();
    const auto& species = system.species();
    const auto& elements = system.elements();

    const unsigned num_phases = phases.size();
    const unsigned bar_size = std::max(unsigned(4), num_phases) * 25;
    const std::string bar1(bar_size, '=');
    const std::string bar2(bar_size, '-');

    unsigned max_size = 0;
    for(const auto& phase : phases)
        max_size = std::max(max_size, phase.numSpecies());

    out << bar1 << std::endl;
    for(const auto& phase : phases)
        out << std::setw(25) << std::left << phase.name();
    out << std::endl;
    out << bar2 << std::endl;
    for(unsigned i = 0; ; ++i)
    {
        if(max_size <= i)
            break;

        for(const auto& phase : phases)
        {
            if(i < phase.numSpecies())
                out << std::setw(25) << std::left << phase.species(i).name();
            else
                out << std::setw(25) << std::left << "";
        }

        out << std::endl;
    }

    out << bar1 << std::endl;
    out << std::setw(25) << std::left << "Index";
    out << std::setw(25) << std::left << "Species";
    out << std::setw(25) << std::left << "Element";
    out << std::setw(25) << std::left << "Phase";
    out << std::endl;
    out << bar2 << std::endl;

    for(unsigned i = 0; ; ++i)
    {
        if(elements.size() <= i && species.size() <= i && phases.size() <= i)
            break;

        out << std::setw(25) << std::left << i;
        out << std::setw(25) << std::left << (i < species.size() ? species[i].name() : "");
        out << std::setw(25) << std::left << (i < elements.size() ? elements[i].name() : "");
        out << std::setw(25) << std::left << (i < phases.size() ? phases[i].name() : "");
        out << std::endl;
    }

    out << bar1 << std::endl;

    return out;
}

} // namespace Reaktoro
