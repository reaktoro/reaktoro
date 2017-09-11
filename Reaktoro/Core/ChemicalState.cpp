// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalProperty.hpp>
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

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

    /// The dual chemical potentials of the elements (in units of J/mol)
    Vector y;

    /// The dual chemical potentials of the species (in units of J/mol)
    Vector z;

    /// Construct a default ChemicalState::Impl instance
    Impl()
    {}

    /// Construct a custom ChemicalState::Impl instance
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialise the molar amounts of the species and dual potentials
        n = zeros(system.numSpecies());
        y = zeros(system.numElements());
        z = zeros(system.numSpecies());
    }

    auto setTemperature(double val) -> void
    {
        Assert(val > 0.0, "Cannot set temperature of the chemical "
            "state with a non-positive value.", "");
        T = val;
    }

    auto setTemperature(double val, std::string units) -> void
    {
        setTemperature(units::convert(val, units, "kelvin"));
    }

    auto setPressure(double val) -> void
    {
        Assert(val > 0.0, "Cannot set pressure of the chemical "
            "state with a non-positive value.", "");
        P = val;
    }

    auto setPressure(double val, std::string units) -> void
    {
        setPressure(units::convert(val, units, "pascal"));
    }

    auto setSpeciesAmounts(double val) -> void
    {
        Assert(val >= 0.0,
            "Cannot set the molar amounts of the species.",
            "The given molar abount is negative.");
        n.fill(val);
    }

    auto setSpeciesAmounts(const Vector& values) -> void
    {
        Assert(static_cast<unsigned>(values.rows()) == system.numSpecies(),
            "Cannot set the molar amounts of the species.",
            "The dimension of the molar abundance vector "
            "is different than the number of species.");
        n = values;
    }

    auto setSpeciesAmounts(const Vector& values, const Indices& indices) -> void
    {
        Assert(static_cast<unsigned>(values.rows()) == indices.size(),
            "Cannot set the molar amounts of the species with given indices.",
            "The dimension of the molar abundance vector "
            "is different than the number of indices.");
        rows(n, indices) = values;
    }

    auto setSpeciesAmount(Index index, double amount) -> void
    {
        Assert(amount >= 0.0,
            "Cannot set the molar amount of the species.",
            "The given molar amount `" + std::to_string(amount) + "` is negative.");
        Assert(index < system.numSpecies(),
            "Cannot set the molar amount of the species.",
            "The given species index is out-of-range.");
        n[index] = amount;
    }

    auto setSpeciesAmount(std::string species, double amount) -> void
    {
        const Index index = system.indexSpeciesWithError(species);
        setSpeciesAmount(index, amount);
    }

    auto setSpeciesAmount(Index index, double amount, std::string units) -> void
    {
        amount = units::convert(amount, units, "mol");
        setSpeciesAmount(index, amount);
    }

    auto setSpeciesAmount(std::string species, double amount, std::string units) -> void
    {
        const Index index = system.indexSpeciesWithError(species);
        setSpeciesAmount(index, amount, units);
    }

    auto setSpeciesMass(Index index, double mass) -> void
    {
        Assert(mass >= 0.0,
            "Cannot set the mass of the species.",
            "The given mass`" + std::to_string(mass) + "` is negative.");
        Assert(index < system.numSpecies(),
            "Cannot set the mass of the species.",
            "The given species index is out-of-range.");
        const double ni = mass/system.species(index).molarMass();
        setSpeciesAmount(index, ni);
    }

    auto setSpeciesMass(std::string species, double mass) -> void
    {
        const Index index = system.indexSpeciesWithError(species);
        setSpeciesMass(index, mass);
    }

    auto setSpeciesMass(Index index, double mass, std::string units) -> void
    {
        mass = units::convert(mass, units, "kg");
        setSpeciesMass(index, mass);
    }

    auto setSpeciesMass(std::string species, double mass, std::string units) -> void
    {
        const Index index = system.indexSpeciesWithError(species);
        setSpeciesMass(index, mass, units);
    }

    /// Set the dual chemical potentials of the species
    auto setSpeciesDualPotentials(const Vector& values) -> void
    {
        Assert(values.size() == z.size(),
            "Could not set the dual chemical potentials of the species.",
            "The dimension of given vector is different from the number of species.");
        z = values;
    }

    /// Set the dual chemical potentials of the elements
    auto setElementDualPotentials(const Vector& values) -> void
    {
        Assert(values.size() == y.size(),
            "Could not set the dual chemical potentials of the elements.",
            "The dimension of given vector is different from the number of elements.");
        y = values;
    }

    auto scaleSpeciesAmounts(double scalar) -> void
    {
        Assert(scalar >= 0.0, "Cannot scale the molar amounts of the species.",
            "The given scalar is negative.");
        for(int i = 0; i < n.rows(); ++i)
            setSpeciesAmount(i, speciesAmount(i) * scalar);
    }

    auto scaleSpeciesAmounts(double scalar, const Indices& indices) -> void
    {
        Assert(scalar >= 0.0, "Cannot scale the molar amounts of the species.",
            "The given scalar is negative.");
        for(Index i : indices)
            setSpeciesAmount(i, speciesAmount(i) * scalar);
    }

    auto scaleSpeciesAmountsInPhase(Index index, double scalar) -> void
    {
        Assert(scalar >= 0.0, "Cannot scale the molar amounts of the species.",
            "The given scalar `" + std::to_string(scalar) << "` is negative.");
        Assert(index < system.numPhases(), "Cannot set the volume of the phase.",
            "The given phase index is out of range.");
        const Index start = system.indexFirstSpeciesInPhase(index);
        const Index size = system.numSpeciesInPhase(index);
        for(unsigned i = 0; i < size; ++i)
            setSpeciesAmount(start + i, speciesAmount(start + i) * scalar);
    }

    auto scalePhaseVolume(Index index, double volume) -> void
    {
        Assert(volume >= 0.0, "Cannot set the volume of the phase.",
            "The given volume is negative.");
        Assert(index < system.numPhases(), "Cannot set the volume of the phase.",
            "The given phase index is out of range.");
        ChemicalProperties properties = system.properties(T, P, n);
        const Vector v = properties.phaseVolumes().val;
        const double scalar = (v[index] != 0.0) ? volume/v[index] : 0.0;
        scaleSpeciesAmountsInPhase(index, scalar);
    }

    auto scalePhaseVolume(Index index, double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scalePhaseVolume(index, volume);
    }

    auto scalePhaseVolume(std::string name, double volume) -> void
    {
        const Index index = system.indexPhase(name);
        scalePhaseVolume(index, volume);
    }

    auto scalePhaseVolume(std::string name, double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scalePhaseVolume(name, volume);
    }

    auto scaleFluidVolume(double volume) -> void
    {
        const auto& fluid_volume = properties().fluidVolume();
        const auto& factor = fluid_volume.val ? volume/fluid_volume.val : 0.0;
        const auto& ifluidspecies = system.indicesFluidSpecies();
        scaleSpeciesAmounts(factor, ifluidspecies);
    }

    auto scaleFluidVolume(double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scaleFluidVolume(volume);
    }

    auto scaleSolidVolume(double volume) -> void
    {
        const auto& solid_volume = properties().solidVolume();
        const auto& factor = solid_volume.val ? volume/solid_volume.val : 0.0;
        const auto& isolidspecies = system.indicesSolidSpecies();
        scaleSpeciesAmounts(factor, isolidspecies);
    }

    auto scaleSolidVolume(double volume, std::string units) -> void
    {
        volume = units::convert(volume, units, "m3");
        scaleSolidVolume(volume);
    }

    auto scaleVolume(double volume) -> void
    {
        Assert(volume >= 0.0, "Cannot set the volume of the chemical state.",
            "The given volume is negative.");
        ChemicalProperties properties = system.properties(T, P, n);
        const Vector v = properties.phaseVolumes().val;
        const double vtotal = sum(v);
        const double scalar = (vtotal != 0.0) ? volume/vtotal : 0.0;
        scaleSpeciesAmounts(scalar);
    }

    auto speciesAmount(Index index) const -> double
    {
        Assert(index < system.numSpecies(),
            "Cannot get the molar amount of the species.",
            "The given index is out-of-range.");
        return n[index];
    }

    auto speciesAmount(std::string name) const -> double
    {
        const Index index = system.indexSpeciesWithError(name);
        return speciesAmount(index);
    }

    auto speciesAmount(Index index, std::string units) const -> double
    {
        return units::convert(speciesAmount(index), "mol", units);
    }

    auto speciesAmount(std::string name, std::string units) const -> double
    {
        const Index index = system.indexSpeciesWithError(name);
        return speciesAmount(index, units);
    }

    auto elementAmounts() const -> Vector
    {
        return system.elementAmounts(n);
    }

    auto elementAmountsInPhase(Index index) const -> Vector
    {
        return system.elementAmountsInPhase(index, n);
    }

    auto elementAmountsInSpecies(const Indices& indices) const -> Vector
    {
        return system.elementAmountsInSpecies(indices, n);
    }

    auto elementAmount(Index ielement) const -> double
    {
        return system.elementAmount(ielement, n);
    }

    auto elementAmount(std::string element) const -> double
    {
        return elementAmount(system.indexElementWithError(element));
    }

    auto elementAmount(Index index, std::string units) const -> double
    {
        return units::convert(elementAmount(index), "mol", units);
    }

    auto elementAmount(std::string name, std::string units) const -> double
    {
        return units::convert(elementAmount(name), "mol", units);
    }

    auto elementAmountInPhase(Index ielement, Index iphase) const -> double
    {
        return system.elementAmountInPhase(ielement, iphase, n);
    }

    auto elementAmountInPhase(std::string element, std::string phase) const -> double
    {
        const unsigned ielement = system.indexElementWithError(element);
        const unsigned iphase = system.indexPhaseWithError(phase);
        return elementAmountInPhase(ielement, iphase);
    }

    auto elementAmountInPhase(Index ielement, Index iphase, std::string units) const -> double
    {
        return units::convert(elementAmountInPhase(ielement, iphase), "mol", units);
    }

    auto elementAmountInPhase(std::string element, std::string phase, std::string units) const -> double
    {
        return units::convert(elementAmountInPhase(element, phase), "mol", units);
    }

    auto elementAmountInSpecies(Index ielement, const Indices& ispecies) const -> double
    {
        return system.elementAmountInSpecies(ielement, ispecies, n);
    }

    auto elementAmountInSpecies(Index ielement, const Indices& ispecies, std::string units) const -> double
    {
        return units::convert(elementAmountInSpecies(ielement, ispecies), "mol", units);
    }

    auto phaseAmount(Index index) const -> double
    {
        const Index first = system.indexFirstSpeciesInPhase(index);
        const Index size = system.numSpeciesInPhase(index);
        return rows(n, first, size).sum();
    }

    auto phaseAmount(std::string name) const -> double
    {
        const Index index = system.indexPhaseWithError(name);
        return phaseAmount(index);
    }

    auto phaseAmount(Index index, std::string units) const -> double
    {
        return units::convert(phaseAmount(index), "mol", units);
    }

    auto phaseAmount(std::string name, std::string units) const -> double
    {
        return units::convert(phaseAmount(name), "mol", units);
    }

    auto properties() const -> ChemicalProperties
    {
        ChemicalProperties res(system);
        res.update(T, P, n);
        return res;
    }

    // Return the stability indices of the phases
    auto phaseStabilityIndices() const -> Vector
    {
        // Auxiliary variables
        const double ln10 = 2.302585092994046;
        const unsigned num_phases = system.numPhases();
        const double RT = universalGasConstant * T;

        // Calculate the normalized z-Lagrange multipliers for all species
        const Vector zRT = z/RT;

        // Initialise the stability indices of the phases
        Vector stability_indices = zeros(num_phases);

        // The index of the first species in each phase iterated below
        unsigned offset = 0;

        // Iterate over all phases
        for(unsigned i = 0 ; i < num_phases; ++i)
        {
            // The number of species in the current phase
            const unsigned num_species = system.numSpeciesInPhase(i);

            if(num_species == 1)
            {
                stability_indices[i] = -zRT[offset]/ln10;
            }
            else
            {
                const Vector zp = rows(zRT, offset, num_species);
                Vector xp = rows(n, offset, num_species);
                const double nsum = sum(xp);
                if(nsum) xp /= nsum; else xp.fill(1.0/num_species);
                stability_indices[i] = std::log10(sum(xp % exp(-zp)));
            }

            offset += num_species;
        }

        return stability_indices;
    }
};

ChemicalState::ChemicalState()
: pimpl(new Impl())
{}

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

auto ChemicalState::setTemperature(double val) -> void
{
    pimpl->setTemperature(val);
}

auto ChemicalState::setTemperature(double val, std::string units) -> void
{
    pimpl->setTemperature(val, units);
}

auto ChemicalState::setPressure(double val) -> void
{
    pimpl->setPressure(val);
}

auto ChemicalState::setPressure(double val, std::string units) -> void
{
    pimpl->setPressure(val, units);
}

auto ChemicalState::setSpeciesAmounts(double val) -> void
{
    pimpl->setSpeciesAmounts(val);
}

auto ChemicalState::setSpeciesAmounts(const Vector& n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::setSpeciesAmounts(const Vector& n, const Indices& indices) -> void
{
    pimpl->setSpeciesAmounts(n, indices);
}

auto ChemicalState::setSpeciesAmount(Index index, double amount) -> void
{
    pimpl->setSpeciesAmount(index, amount);
}

auto ChemicalState::setSpeciesAmount(std::string species, double amount) -> void
{
    pimpl->setSpeciesAmount(species, amount);
}

auto ChemicalState::setSpeciesAmount(Index index, double amount, std::string units) -> void
{
    pimpl->setSpeciesAmount(index, amount, units);
}

auto ChemicalState::setSpeciesAmount(std::string species, double amount, std::string units) -> void
{
    pimpl->setSpeciesAmount(species, amount, units);
}

auto ChemicalState::setSpeciesMass(Index index, double mass) -> void
{
    pimpl->setSpeciesMass(index, mass);
}

auto ChemicalState::setSpeciesMass(std::string name, double mass) -> void
{
    pimpl->setSpeciesMass(name, mass);
}

auto ChemicalState::setSpeciesMass(Index index, double mass, std::string units) -> void
{
    pimpl->setSpeciesMass(index, mass, units);
}

auto ChemicalState::setSpeciesMass(std::string name, double mass, std::string units) -> void
{
    pimpl->setSpeciesMass(name, mass, units);
}

auto ChemicalState::setSpeciesDualPotentials(const Vector& z) -> void
{
    pimpl->setSpeciesDualPotentials(z);
}

auto ChemicalState::setElementDualPotentials(const Vector& y) -> void
{
    pimpl->setElementDualPotentials(y);
}

auto ChemicalState::scaleSpeciesAmounts(double scalar) -> void
{
    pimpl->scaleSpeciesAmounts(scalar);
}

auto ChemicalState::scaleSpeciesAmountsInPhase(Index index, double scalar) -> void
{
    pimpl->scaleSpeciesAmountsInPhase(index, scalar);
}

auto ChemicalState::scalePhaseVolume(Index index, double volume) -> void
{
    pimpl->scalePhaseVolume(index, volume);
}

auto ChemicalState::scalePhaseVolume(Index index, double volume, std::string units) -> void
{
    pimpl->scalePhaseVolume(index, volume, units);
}

auto ChemicalState::scalePhaseVolume(std::string name, double volume) -> void
{
    pimpl->scalePhaseVolume(name, volume);
}

auto ChemicalState::scalePhaseVolume(std::string name, double volume, std::string units) -> void
{
    pimpl->scalePhaseVolume(name, volume, units);
}

auto ChemicalState::scaleFluidVolume(double volume) -> void
{
    pimpl->scaleFluidVolume(volume);
}

auto ChemicalState::scaleFluidVolume(double volume, std::string units) -> void
{
    pimpl->scaleFluidVolume(volume, units);
}

auto ChemicalState::scaleSolidVolume(double volume) -> void
{
    pimpl->scaleSolidVolume(volume);
}

auto ChemicalState::scaleSolidVolume(double volume, std::string units) -> void
{
    pimpl->scaleSolidVolume(volume, units);
}

auto ChemicalState::scaleVolume(double volume) -> void
{
    pimpl->scaleVolume(volume);
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

auto ChemicalState::speciesAmounts(const Indices& indices) const -> Vector
{
    return rows(pimpl->n, indices);
}

auto ChemicalState::speciesAmount(Index index) const -> double
{
    return pimpl->speciesAmount(index);
}

auto ChemicalState::speciesAmount(std::string name) const -> double
{
    return pimpl->speciesAmount(name);
}

auto ChemicalState::speciesAmount(Index index, std::string units) const -> double
{
    return pimpl->speciesAmount(index, units);
}

auto ChemicalState::speciesAmount(std::string species, std::string units) const -> double
{
    return pimpl->speciesAmount(species, units);
}

auto ChemicalState::speciesDualPotentials() const -> const Vector&
{
    return pimpl->z;
}

auto ChemicalState::elementAmounts() const -> Vector
{
    return pimpl->elementAmounts();
}

auto ChemicalState::elementAmountsInPhase(Index iphase) const -> Vector
{
    return pimpl->elementAmountsInPhase(iphase);
}

auto ChemicalState::elementAmountsInSpecies(const Indices& ispecies) const -> Vector
{
    return pimpl->elementAmountsInSpecies(ispecies);
}

auto ChemicalState::elementAmount(Index ielement) const -> double
{
    return pimpl->elementAmount(ielement);
}

auto ChemicalState::elementAmount(std::string element) const -> double
{
    return pimpl->elementAmount(element);
}

auto ChemicalState::elementAmount(Index index, std::string units) const -> double
{
    return pimpl->elementAmount(index, units);
}

auto ChemicalState::elementAmount(std::string name, std::string units) const -> double
{
    return pimpl->elementAmount(name, units);
}

auto ChemicalState::elementAmountInPhase(Index ielement, Index iphase) const -> double
{
    return pimpl->elementAmountInPhase(ielement, iphase);
}

auto ChemicalState::elementAmountInPhase(std::string element, std::string phase) const -> double
{
    return pimpl->elementAmountInPhase(element, phase);
}

auto ChemicalState::elementAmountInPhase(Index ielement, Index iphase, std::string units) const -> double
{
    return pimpl->elementAmountInPhase(ielement, iphase, units);
}

auto ChemicalState::elementAmountInPhase(std::string element, std::string phase, std::string units) const -> double
{
    return pimpl->elementAmountInPhase(element, phase, units);
}

auto ChemicalState::elementAmountInSpecies(Index ielement, const Indices& ispecies) const -> double
{
    return pimpl->elementAmountInSpecies(ielement, ispecies);
}

auto ChemicalState::elementAmountInSpecies(Index ielement, const Indices& ispecies, std::string units) const -> double
{
    return pimpl->elementAmountInSpecies(ielement, ispecies, units);
}

auto ChemicalState::elementDualPotentials() const -> const Vector&
{
    return pimpl->y;
}

auto ChemicalState::phaseAmount(Index index) const -> double
{
    return pimpl->phaseAmount(index);
}

auto ChemicalState::phaseAmount(std::string name) const -> double
{
    return pimpl->phaseAmount(name);
}

auto ChemicalState::phaseAmount(Index index, std::string units) const -> double
{
    return pimpl->phaseAmount(index, units);
}

auto ChemicalState::phaseAmount(std::string name, std::string units) const -> double
{
    return pimpl->phaseAmount(name, units);
}

auto ChemicalState::phaseStabilityIndices() const -> Vector
{
    return pimpl->phaseStabilityIndices();
}

auto ChemicalState::properties() const -> ChemicalProperties
{
    return pimpl->properties();
}

auto ChemicalState::output(std::string filename) -> void
{
    std::ofstream out(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&
{
    const ChemicalSystem& system = state.system();
    const double& T = state.temperature();
    const double& P = state.pressure();
    const double& R = universalGasConstant;
    const double& F = faradayConstant;
    const Vector& n = state.speciesAmounts();
    const Vector& y = state.elementDualPotentials();
    const Vector& z = state.speciesDualPotentials();
    const ChemicalProperties properties = state.properties();
    const Vector molar_fractions = properties.moleFractions().val;
    const Vector activity_coeffs = exp(properties.lnActivityCoefficients().val);
    const Vector activities = exp(properties.lnActivities().val);
    const Vector chemical_potentials = properties.chemicalPotentials().val;
    const Vector phase_moles = properties.phaseAmounts().val;
    const Vector phase_masses = properties.phaseMasses().val;
    const Vector phase_molar_volumes = properties.phaseMolarVolumes().val;
    const Vector phase_volumes = properties.phaseVolumes().val;
    const Vector phase_volume_fractions = phase_volumes/sum(phase_volumes);
    const Vector phase_densities = phase_masses/phase_volumes;
    const Vector phase_stability_indices = state.phaseStabilityIndices();

    // Calculate pH, pE, and Eh
    const double I  = ChemicalProperty::ionicStrength(system)(properties).val;
    const double pH = ChemicalProperty::pH(system)(properties).val;
    const double pE = ChemicalProperty::pE(system)(properties).val;
    const double Eh = std::log(10)*R*T/F*pE;
    const double alk = ChemicalProperty::alkalinity(system)(properties).val;

    const unsigned num_phases = system.numPhases();
    const unsigned bar_size = std::max(unsigned(9), num_phases + 2) * 25;
    const std::string bar1(bar_size, '=');
    const std::string bar2(bar_size, '-');

    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Temperature [K]";
    out << std::left << std::setw(25) << "Temperature [°C]";
    out << std::left << std::setw(25) << "Pressure [MPa]";
    out << std::endl << bar2 << std::endl;

    out << std::left << std::setw(25) << T;
    out << std::left << std::setw(25) << T - 273.15;
    out << std::left << std::setw(25) << P * 1e-6;
    out << std::endl;

    // Set output in scientific notation
    auto flags = out.flags();
    out << std::setprecision(6);

    // Output the table of the element-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Element";
    out << std::left << std::setw(25) << "Amount [mol]";
    for(const auto& phase : system.phases())
        out << std::left << std::setw(25) << phase.name() + " [mol]";
    out << std::left << std::setw(25) << "Dual Potential [kJ/mol]";
    out << std::endl;
    out << bar2 << std::endl;
    for(unsigned i = 0; i < system.numElements(); ++i)
    {
        out << std::left << std::setw(25) << system.element(i).name();
        out << std::left << std::setw(25) << system.elementAmount(i, n);
        for(unsigned j = 0; j < system.numPhases(); ++j)
            out << std::left << std::setw(25) << system.elementAmountInPhase(i, j, n);
        out << std::left << std::setw(25) << y[i]/1000; // convert from J/mol to kJ/mol
        out << std::endl;
    }

    // Output the table of the species-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Species";
    out << std::left << std::setw(25) << "Amount [mol]";
    out << std::left << std::setw(25) << "Mole Fraction [mol/mol]";
    out << std::left << std::setw(25) << "Activity Coefficient [-]";
    out << std::left << std::setw(25) << "Activity [-]";
    out << std::left << std::setw(25) << "Potential [kJ/mol]";
    out << std::left << std::setw(25) << "Dual Potential [kJ/mol]";
    out << std::endl;
    out << bar2 << std::endl;
    for(unsigned i = 0; i < system.numSpecies(); ++i)
    {
        out << std::left << std::setw(25) << system.species(i).name();
        out << std::left << std::setw(25) << n[i];
        out << std::left << std::setw(25) << molar_fractions[i];
        out << std::left << std::setw(25) << activity_coeffs[i];
        out << std::left << std::setw(25) << activities[i];
        out << std::left << std::setw(25) << chemical_potentials[i]/1000; // convert from J/mol to kJ/mol
        out << std::left << std::setw(25) << z[i]/1000; // convert from J/mol to kJ/mol
        out << std::endl;
    }

    // Output the table of the phase-related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Phase";
    out << std::left << std::setw(25) << "Amount [mol]";
    out << std::left << std::setw(25) << "Stability";
    out << std::left << std::setw(25) << "Stability Index [-]";
    out << std::left << std::setw(25) << "Mass [kg]";
    out << std::left << std::setw(25) << "Volume [m3]";
    out << std::left << std::setw(25) << "Density [kg/m3]";
    out << std::left << std::setw(25) << "Molar Volume [m3/mol]";
    out << std::left << std::setw(25) << "Volume Fraction [m3/m3]";
    out << std::endl;
    out << bar2 << std::endl;
    for(unsigned i = 0; i < system.numPhases(); ++i)
    {
        int extra = (phase_stability_indices[i] < 0 ? 0 : 1);
        std::string stability = std::abs(phase_stability_indices[i]) < 1e-2 ? "stable" : "unstable";
        out << std::left << std::setw(25) << system.phase(i).name();
        out << std::left << std::setw(25) << phase_moles[i];
        out << std::setw(25 + extra) << std::left << stability;
        out << std::setw(25 - extra) << std::left << phase_stability_indices[i];
        out << std::left << std::setw(25) << phase_masses[i];
        out << std::left << std::setw(25) << phase_volumes[i];
        out << std::left << std::setw(25) << phase_densities[i];
        out << std::left << std::setw(25) << phase_molar_volumes[i];
        out << std::left << std::setw(25) << phase_volume_fractions[i];
        out << std::endl;
    }

    // Output the table of the aqueous phase related state
    out << bar1 << std::endl;
    out << std::left << std::setw(25) << "Ionic Strength [molal]";
    out << std::left << std::setw(25) << "pH";
    out << std::left << std::setw(25) << "pE";
    out << std::left << std::setw(25) << "Reduction Potential [V]";
    out << std::left << std::setw(25) << "Alkalinity [eq/L]";
    out << std::endl << bar2 << std::endl;
    out << std::left << std::setw(25) << I;
    out << std::left << std::setw(25) << pH;
    out << std::left << std::setw(25) << pE;
    out << std::left << std::setw(25) << Eh;
    out << std::left << std::setw(25) << alk;
    out << std::endl << bar1 << std::endl;

    // Recover the previous state of `out`
    out.flags(flags);

    return out;
}

auto operator+(const ChemicalState& l, const ChemicalState& r) -> ChemicalState
{
    const Vector& nl = l.speciesAmounts();
    const Vector& nr = r.speciesAmounts();
    ChemicalState res = l;
    res.setSpeciesAmounts(nl + nr);
    return res;
}

auto operator*(double scalar, const ChemicalState& state) -> ChemicalState
{
    ChemicalState res = state;
    res.scaleSpeciesAmounts(scalar);
    return res;
}

auto operator*(const ChemicalState& state, double scalar) -> ChemicalState
{
    return scalar*state;
}

} // namespace Reaktoro
