/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ChemicalState.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Utils/ElementUtils.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {
namespace internal {

inline auto checkSpeciesDimension(const ChemicalSystem& system, const Vector& n) -> void
{
    if(system.numSpecies() != n.rows())
    error("Cannot set the composition of the species with the given vector of species moles.",
        "The dimension of the provided vector does not match the number of species.");
}

inline auto checkEquilibriumSpeciesDimension(const Partitioning& partitioning, const Vector& ne) -> void
{
    if(partitioning.numEquilibriumSpecies() != ne.rows())
    error("Cannot set the composition of the equilibrium species with the given vector of species moles.",
        "The dimension of the provided vector does not match the number of equilibrium species.");
}

inline auto checkKineticSpeciesDimension(const Partitioning& partitioning, const Vector& nk) -> void
{
    if(partitioning.numKineticSpecies() != nk.rows())
    error("Cannot set the composition of the kinetic species with the given vector of species moles.",
        "The dimension of the provided vector does not match the number of kinetic species.");
}

inline auto checkInertSpeciesDimension(const Partitioning& partitioning, const Vector& ni) -> void
{
    if(partitioning.numInertSpecies() != ni.rows())
    error("Cannot set the composition of the inert species with the given vector of species moles.",
        "The dimension of the provided vector does not match the number of inert species.");
}

inline auto checkPorosityValue(double porosity) -> void
{
    if(porosity < 0.0 or porosity > 1.0)
    error("Cannot set the composition of the minerals with given porosity value.",
        "The value provided is not between 0 and 1.");
}

inline auto nonAmountOrMassUnitError(const std::string& unit) -> void
{
    Exception exception;
    exception.error << "Cannot set the amount or mass of a species.";
    exception.reason << "The provided unit " << unit << " is not convertible to units of amount or mass (e.g., mol or g).";
    raise(exception);
}

inline auto mineralVolumePercentError(const std::string& composition) -> void
{
    Exception exception;
    exception.error << "Cannot set the composition of the minerals.";
    exception.reason << "The volume percent of the minerals in the composition string `" << composition << " does not sum to 100.";
    raise(exception);
}

} /* namespace internal */

using namespace internal;

class ChemicalState::Impl
{
public:
    /// The chemical system instance
    ChemicalSystem system;

    /// The temperature state of the chemical system
    double T = 298.15;

    /// The pressure state of the chemical system
    double P = 1.0e+05;

    /// The number of moles of every species in the chemical system
    Vector n;

public:
    explicit Impl(const ChemicalSystem& system)
    : system(system), n(system.numSpecies())
    {
        setTemperature(25 * unit(degC));
        setPressure(1 * unit(bar));
        setComposition(1.0e-7);
    }

    auto setTemperature(units::Temperature value) -> void
    {
        T = value.in(unit(K));
    }

    auto setTemperature(double value, const std::string& unit) -> void
    {
        setTemperature(units::Temperature(value, unit));
    }

    auto setPressure(units::Pressure value) -> void
    {
        P = value.in(unit(Pa));
    }

    auto setPressure(double value, const std::string& unit) -> void
    {
        setPressure(units::Pressure(value, unit));
    }

    auto setComposition(double value) -> void
    {
        n.fill(value);
    }

    auto setComposition(const Vector& n) -> void
    {
        checkSpeciesDimension(system, n);
        this->n = n;
    }

    auto setComposition(const Indices& indices, const Vector& values) -> void
    {
        setRows(indices, values, n);
    }

    auto setEquilibriumComposition(const Partitioning& partitioning, const Vector& ne) -> void
    {
        checkEquilibriumSpeciesDimension(partitioning, ne);
        partitioning.setEquilibriumRows(ne, n);
    }

    auto setKineticComposition(const Partitioning& partitioning, const Vector& nk) -> void
    {
        checkKineticSpeciesDimension(partitioning, nk);
        partitioning.setKineticRows(nk, n);
    }

    auto setInertComposition(const Partitioning& partitioning, const Vector& ni) -> void
    {
        checkInertSpeciesDimension(partitioning, ni);
        partitioning.setInertRows(ni, n);
    }

    auto setSpeciesAmount(Index ispecies, units::Amount value) -> void
    {
        n[ispecies] = value.in(unit(mol));
    }

    auto setSpeciesAmount(const std::string& species, units::Amount value) -> void
    {
        auto ispecies = system.idxSpeciesWithError(species);
        setSpeciesAmount(ispecies, value);
    }

    auto setSpeciesMass(Index ispecies, units::Mass value) -> void
    {
        const units::MolarMass molar_mass = system.species(ispecies).molarMass();
        const units::Amount moles = value/molar_mass;
        setSpeciesAmount(ispecies, moles);
    }

    auto setSpeciesMass(const std::string& species, units::Mass value) -> void
    {
        auto ispecies = system.idxSpeciesWithError(species);
        setSpeciesMass(ispecies, value);
    }

    auto setMinerals(const std::string& composition, double volume, const std::string& unit, double porosity) -> void
    {
        checkPorosityValue(porosity);

        // The total volume of the minerals (in units of m3)
        const double vol = units::convert(volume, unit, "m3") * (1.0 - porosity);

        // The auxiliary variable to check later if the sum of each volume percent equals 100
        double total = 0.0;

        // Iterate over all pairs "xm:mineral" in the composition string "x1:mineral1 x2:mineral2 ..."
        for(auto word : split(composition, " "))
        {
            // Get the pair (x, mineral) from the current word "x:mineral"
            const auto pair    = split(word, ":");
            const auto xm      = 0.01 * tofloat(pair.front());
            const auto mineral = pair.back();

            // Update the `total` variable
            total += xm;

            // Calculate the molar volume of the mineral (in units of m3/mol)
            const auto vm = system.species(mineral).molarVolume(T, P);

            // The number of moles of the mineral
            const auto nm = xm*vol/vm;

            // Set the amount of the mineral
            setSpeciesAmount(mineral, nm * unit(mol));
        }

        // Check if the total of the volume percent of each mineral sums to 100
        if(std::abs(total - 1.0) > 1e-6)
            mineralVolumePercentError(composition);
    }

    auto set(const std::string& species, double value, const std::string& unit) -> void
    {
        if(units::convertible(unit, "mol"))
            setSpeciesAmount(species, units::Amount(value, unit));
        else if(units::convertible(unit, "g"))
            setSpeciesMass(species, units::Mass(value, unit));
        else nonAmountOrMassUnitError(unit);
    }

    auto phaseComposition(Index iphase) const -> Vector
    {
        return rows(system.idxSpeciesInPhase(iphase), n);
    }

    auto equilibriumComposition(const Partitioning& partitioning) -> Vector
    {
        return partitioning.equilibriumRows(n);
    }

    auto kineticComposition(const Partitioning& partitioning) -> Vector
    {
        return partitioning.kineticRows(n);
    }

    auto inertComposition(const Partitioning& partitioning) -> Vector
    {
        return partitioning.inertRows(n);
    }

    auto amountElement(const Index& idx_element) const -> units::Amount
    {
        const std::string element = system.element(idx_element);

        return amountElement(element);
    }

    auto amountElement(const std::string& element) const -> units::Amount
    {
        double total = 0.0;
        for(Index i = 0; i < system.numSpecies(); ++i)
            total += n[i] * system.species(i).elementAtoms(element);
        return total * unit(mol);
    }

    auto amountElement(const std::string& element, const std::string& phase) const -> units::Amount
    {
        double total = 0.0;
        const auto iphase = system.idxPhase(phase);
        for(Index i : system.idxSpeciesInPhase(iphase))
            total += n[i] * system.species(i).elementAtoms(element);
        return total * unit(mol);
    }

    auto amountPhase(const Index& idx_phase) const -> units::Amount
    {
        return rows(system.idxSpeciesInPhase(idx_phase), n).sum() * unit(mol);
    }

    auto amountPhase(const std::string& phase) const -> units::Amount
    {
        return amountPhase(system.idxPhase(phase));
    }

    auto massElement(const std::string& element) const -> units::Mass
    {
        return atomicMass(element) * amountElement(element);
    }

    auto massElement(const std::string& element, const std::string& phase) const -> units::Mass
    {
        return atomicMass(element) * amountElement(element, phase);
    }

    auto molalityElement(const std::string& element) const -> units::Molality
    {
        const double ne = amountElement(element, "Aqueous").in(unit(mol));;
        const double nw = amount("H2O(l)").in(unit(mol));
        const double mw = mass("H2O(l)").in(unit(kg));
        const double value = (element == "H") ? (ne - 2*nw)/mw : (element == "O") ? (ne - nw)/mw : ne/mw;
        return value * unit(molal);
    }

    auto molality(const std::string& species) const -> units::Molality
    {
        return amount(species)/mass("H2O(l)");
    }

    auto amount(Index ispecies) const -> units::Amount
    {
        const double value = ispecies > system.numSpecies() ? 0.0 : n[ispecies];
        return value * unit(mol);
    }

    auto amount(const std::string& species) const -> units::Amount
    {
        const auto ispecies = system.idxSpeciesWithError(species);
        return amount(ispecies);
    }

    auto mass(Index ispecies) const -> units::Mass
    {
        const auto molar_mass = system.species(ispecies).molarMass();
        return amount(ispecies) * molar_mass;
    }

    auto mass(const std::string& species) const -> units::Mass
    {
        const auto ispecies = system.idxSpeciesWithError(species);
        return mass(ispecies);
    }

    auto molarFractions() const -> Vector
    {
        return system.molarFractions(n);
    }

    auto concentrations() const -> Vector
    {
        return system.concentrations(n);
    }

    auto activities() const -> PartialVector
    {
        return system.activities(T, P, n);
    }

    auto activity(const PartialVector& a, Index ispecies) const -> double
    {
        return func(a)[ispecies];
    }

    auto activity(const PartialVector& a, const std::string& species) const -> double
    {
        const auto ispecies = system.idxSpeciesWithError(species);
        return activity(a, ispecies);
    }

    auto acidity(const PartialVector& a) const -> double
    {
        return -std::log10(activity(a, "H+"));
    }
};

ChemicalState::ChemicalState(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalState::ChemicalState(const ChemicalState& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalState::~ChemicalState()
{}

auto ChemicalState::operator=(ChemicalState other) -> ChemicalState&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalState::temperature() const -> double
{
    return pimpl->T;
}

auto ChemicalState::pressure() const -> double
{
    return pimpl->P;
}

auto ChemicalState::composition() const -> const Vector&
{
    return pimpl->n;
}

auto ChemicalState::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalState::setTemperature(units::Temperature value) -> void
{
    pimpl->setTemperature(value);
}

auto ChemicalState::setTemperature(double value, const std::string& unit) -> void
{
    pimpl->setTemperature(value, unit);
}

auto ChemicalState::setPressure(units::Pressure value) -> void
{
    pimpl->setPressure(value);
}

auto ChemicalState::setPressure(double value, const std::string& unit) -> void
{
    pimpl->setPressure(value, unit);
}

auto ChemicalState::setComposition(double value) -> void
{
    pimpl->setComposition(value);
}

auto ChemicalState::setComposition(const Vector& n) -> void
{
    pimpl->setComposition(n);
}

auto ChemicalState::setComposition(const Indices& indices, const Vector& values) -> void
{
    pimpl->setComposition(indices, values);
}

auto ChemicalState::setEquilibriumComposition(const Partitioning& partitioning, const Vector& ne) -> void
{
    pimpl->setEquilibriumComposition(partitioning, ne);
}

auto ChemicalState::setKineticComposition(const Partitioning& partitioning, const Vector& nk) -> void
{
    pimpl->setKineticComposition(partitioning, nk);
}

auto ChemicalState::setInertComposition(const Partitioning& partitioning, const Vector& ni) -> void
{
    pimpl->setInertComposition(partitioning, ni);
}

auto ChemicalState::setSpeciesAmount(Index ispecies, units::Amount value) -> void
{
    pimpl->setSpeciesAmount(ispecies, value);
}

auto ChemicalState::setSpeciesAmount(const std::string& species, units::Amount value) -> void
{
    pimpl->setSpeciesAmount(species, value);
}

auto ChemicalState::setSpeciesMass(Index ispecies, units::Mass value) -> void
{
    pimpl->setSpeciesMass(ispecies, value);
}

auto ChemicalState::setSpeciesMass(const std::string& species, units::Mass value) -> void
{
    pimpl->setSpeciesMass(species, value);
}

auto ChemicalState::setMinerals(const std::string& composition, double volume, const std::string& unit, double porosity) -> void
{
    pimpl->setMinerals(composition, volume, unit, porosity);
}

auto ChemicalState::set(const std::string& species, double value, const std::string& unit) -> void
{
    pimpl->set(species, value, unit);
}

auto ChemicalState::phaseComposition(Index iphase) const -> Vector
{
    return pimpl->phaseComposition(iphase);
}

auto ChemicalState::equilibriumComposition(const Partitioning& partitioning) const -> Vector
{
    return pimpl->equilibriumComposition(partitioning);
}

auto ChemicalState::kineticComposition(const Partitioning& partitioning) const -> Vector
{
    return pimpl->kineticComposition(partitioning);
}

auto ChemicalState::inertComposition(const Partitioning& partitioning) const -> Vector
{
    return pimpl->inertComposition(partitioning);
}

auto ChemicalState::amount(Index ispecies) const -> units::Amount
{
    return pimpl->amount(ispecies);
}

auto ChemicalState::amount(const std::string& species) const -> units::Amount
{
    return pimpl->amount(species);
}

auto ChemicalState::amountElement(const Index& idx_element) const -> units::Amount
{
    return pimpl->amountElement(idx_element);
}

auto ChemicalState::amountElement(const std::string& element) const -> units::Amount
{
    return pimpl->amountElement(element);
}

auto ChemicalState::amountElement(const std::string& element, const std::string& phase) const -> units::Amount
{
    return pimpl->amountElement(element, phase);
}

auto ChemicalState::amountPhase(const Index& idx_phase) const -> units::Amount
{
    return pimpl->amountPhase(idx_phase);
}

auto ChemicalState::amountPhase(const std::string& phase) const -> units::Amount
{
    return pimpl->amountPhase(phase);
}

auto ChemicalState::mass(Index ispecies) const -> units::Mass
{
    return pimpl->mass(ispecies);
}

auto ChemicalState::mass(const std::string& species) const -> units::Mass
{
    return pimpl->mass(species);
}

auto ChemicalState::massElement(const std::string& element) const -> units::Mass
{
    return pimpl->massElement(element);
}

auto ChemicalState::massElement(const std::string& element, const std::string& phase) const -> units::Mass
{
    return pimpl->massElement(element, phase);
}

auto ChemicalState::molality(const std::string& species) const -> units::Molality
{
    return pimpl->molality(species);
}

auto ChemicalState::molalityElement(const std::string& element) const -> units::Molality
{
    return pimpl->molalityElement(element);
}

auto ChemicalState::molarFractions() const -> Vector
{
    return pimpl->molarFractions();
}

auto ChemicalState::concentrations() const -> Vector
{
    return pimpl->concentrations();
}

auto ChemicalState::activities() const -> PartialVector
{
    return pimpl->activities();
}

auto ChemicalState::activity(const PartialVector& a, Index ispecies) const -> double
{
    return pimpl->activity(a, ispecies);
}

auto ChemicalState::activity(const PartialVector& a, const std::string& species) const -> double
{
    return pimpl->activity(a, species);
}

auto ChemicalState::acidity(const PartialVector& a) const -> double
{
    return pimpl->acidity(a);
}

auto ChemicalState::acidity() const -> double
{
    return acidity(activities());
}

auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&
{
    // Auxiliary reference to the chemical system
    const ChemicalSystem& system = state.system();

    // Auxiliary variables for pretty-printing
    const unsigned nfill = 5 * 25;
    const std::string bar(nfill, '=');

    // The molar composition, activity, activity coefficients, and concentrations of the species
    Vector n = state.composition();
    Vector a = func(state.activities());
    Vector c = state.concentrations();
    Vector g = a.array() / c.array();

    // Fix the activity and activity coefficients of species that have zero number of moles
    a = (n.array() > 0.0).select(a, 0.0);
    g = (n.array() > 0.0).select(g, 1.0);

    // Correct the activity coefficient value of the water species
    const Index iH2O = system.idxSpeciesWithError("H2O(l)");
    if(iH2O < system.numSpecies())
        g[iH2O] = a[iH2O];

    // Output the header
    out << bar << std::endl;
    out << std::setw(25) << std::left << "Species";
    out << std::setw(25) << std::left << "Amount";
    out << std::setw(25) << std::left << "Activity";
    out << std::setw(25) << std::left << "Activity Coefficient";
    out << std::setw(25) << std::left << "Concentrations" << std::endl;
    out << bar << std::endl;

    // Output the chemical state
    for(unsigned i = 0; i < state.system().numSpecies(); ++i)
    {
        const std::string name = state.system().species(i).name();

        out << std::setw(25) << std::left << name;
        out << std::setw(25) << std::left << n[i];
        out << std::setw(25) << std::left << a[i];
        out << std::setw(25) << std::left << g[i];
        out << std::setw(25) << std::left << c[i];
        out << std::endl;
    }

    return out;
}

} // namespace Reaktor
