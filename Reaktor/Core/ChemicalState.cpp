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

#include "ChemicalState.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {

struct ChemicalState::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The temperature state of the chemical system (in units of K)
    double T = 298.15;

    /// The pressure state of the chemical system (in units of Pa)
    double P = 1.0e+05;

    /// The molar amounts of the chemical species
    Vector n;

    /// The Lagrange multiplier with respect to the equilibrium charge balance constraint
    double ycharge = 0.0;

    /// The Lagrange multipliers with respect to the equilibrium mass balance constraints (in units of J/mol)
    Vector y;

    /// The Lagrange multipliers with respect to the bound constraints of the species (in units of J/mol)
    Vector z;

    /// Construct a custom ChemicalState::Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        n = zeros(system.numSpecies());
        y = zeros(system.numElements());
        z = zeros(system.numSpecies());
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

auto ChemicalState::setTemperature(double val) const -> void
{
    Assert(val > 0.0, "Cannot set temperature of the chemical state with a non-positive value.", "");
    pimpl->T = val;
}

auto ChemicalState::setTemperature(double val, std::string units) const -> void
{
    setTemperature(units::convert(val, units, "kelvin"));
}

auto ChemicalState::setPressure(double val) const -> void
{
    Assert(val > 0.0, "Cannot set pressure of the chemical state with a non-positive value.", "");
    pimpl->P = val;
}

auto ChemicalState::setPressure(double val, std::string units) const -> void
{
    setPressure(units::convert(val, units, "pascal"));
}

auto ChemicalState::setSpeciesAmounts(double val) -> void
{
    Assert(val >= 0.0,
        "Cannot set the molar amounts of the species.",
        "The given molar abount is negative.");
    pimpl->n.fill(val);
}

auto ChemicalState::setSpeciesAmounts(const Vector& n) -> void
{
    Assert(n.rows() == system().numSpecies(),
        "Cannot set the molar amounts of the species.",
        "The dimension of the molar abundance vector "
        "is different than the number of species.");
    pimpl->n = n;
}

auto ChemicalState::setSpeciesAmount(Index index, double amount) -> void
{
    const unsigned num_species = system().species().size();
    Assert(amount >= 0.0,
        "Cannot set the molar amount of the species.",
        "The given molar amount is negative.");
    Assert(index < num_species,
        "Cannot set the molar amount of the species.",
        "The given index is out-of-range.");
    pimpl->n[index] = amount;
}

auto ChemicalState::setSpeciesAmount(std::string species, double amount) -> void
{
    const Index index = system().indexSpeciesWithError(species);
    setSpeciesAmount(index, amount);
}

auto ChemicalState::setSpeciesAmount(Index index, double amount, std::string units) -> void
{
    setSpeciesAmount(index, units::convert(amount, units, "mol"));
}

auto ChemicalState::setSpeciesAmount(std::string species, double amount, std::string units) -> void
{
    setSpeciesAmount(species, units::convert(amount, units, "mol"));
}

auto ChemicalState::setChargePotential(double ycharge) -> void
{
    pimpl->ycharge = ycharge;
}

auto ChemicalState::setElementPotentials(const Vector& y) -> void
{
    pimpl->y = y;
}

auto ChemicalState::setSpeciesPotentials(const Vector& z) -> void
{
    pimpl->z = z;
}

auto ChemicalState::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalState::temperature() const -> double
{
    return pimpl->T;
}

auto ChemicalState::pressure() const -> double
{
    return pimpl->P;
}

auto ChemicalState::speciesAmounts() const -> const Vector&
{
    return pimpl->n;
}

auto ChemicalState::chargePotential() const -> double
{
    return pimpl->ycharge;
}

auto ChemicalState::elementPotentials() const -> const Vector&
{
    return pimpl->y;
}

auto ChemicalState::speciesPotentials() const -> const Vector&
{
    return pimpl->z;
}

auto ChemicalState::speciesAmount(Index index) const -> double
{
    Assert(index < system().numSpecies(),
        "Cannot get the molar amount of the species.",
        "The given index is out-of-range.");
    return pimpl->n[index];
}

auto ChemicalState::speciesAmount(std::string name) const -> double
{
    return speciesAmount(system().indexSpeciesWithError(name));
}

auto ChemicalState::speciesAmount(Index ispecies, std::string units) const -> double
{
    return units::convert(speciesAmount(ispecies), "mol", units);
}

auto ChemicalState::speciesAmount(std::string species, std::string units) const -> double
{
    return units::convert(speciesAmount(species), "mol", units);
}

auto ChemicalState::speciesAmountsInPhase(Index index) const -> Vector
{
    return system().nphase(pimpl->n, index);
}

auto ChemicalState::speciesAmountsInPhase(std::string name) const -> Vector
{
    return speciesAmountsInPhase(system().indexPhaseWithError(name));
}

auto ChemicalState::elementAmounts() const -> Vector
{
    const Matrix& W = system().formulaMatrix();
    return W * pimpl->n;
}

auto ChemicalState::elementAmount(Index ielement) const -> double
{
    const Matrix& W = system().formulaMatrix();
    return W.row(ielement) * pimpl->n;
}

auto ChemicalState::elementAmount(std::string element) const -> double
{
    return elementAmount(system().indexElementWithError(element));
}

auto ChemicalState::elementAmount(Index index, std::string units) const -> double
{
    return units::convert(elementAmount(index), "mol", units);
}

auto ChemicalState::elementAmount(std::string name, std::string units) const -> double
{
    return units::convert(elementAmount(name), "mol", units);
}

auto ChemicalState::elementAmountInPhase(Index ielement, Index iphase) const -> double
{
    const Matrix& W = system().formulaMatrix();
    const unsigned offset = system().offset(iphase);
    const unsigned size = system().phase(iphase).numSpecies();
    const auto Wp = cols(W, offset, size);
    const auto np = cols(pimpl->n, offset, size);
    return dot(Wp.row(ielement), np);
}

auto ChemicalState::elementAmountInPhase(std::string element, std::string phase) const -> double
{
    const unsigned ielement = system().indexElementWithError(element);
    const unsigned iphase = system().indexPhaseWithError(phase);
    return elementAmountInPhase(ielement, iphase);
}

auto ChemicalState::elementAmountInPhase(Index ielement, Index iphase, std::string units) const -> double
{
    return units::convert(elementAmountInPhase(ielement, iphase), "mol", units);
}

auto ChemicalState::elementAmountInPhase(std::string element, std::string phase, std::string units) const -> double
{
    return units::convert(elementAmountInPhase(element, phase), "mol", units);
}

auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&
{
    return out;
}

} // namespace Reaktor
