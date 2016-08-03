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

#include "EquilibriumState.hpp"

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>

namespace Reaktoro {

struct EquilibriumState::Impl
{
    /// The dual chemical potentials of the elements (in units of J/mol)
    Vector y;

    /// The dual chemical potentials of the species (in units of J/mol)
    Vector z;

    /// Construct a default EquilibriumState::Impl instance
    Impl()
    {}

    /// Construct a custom EquilibriumState::Impl instance
    Impl(const ChemicalSystem& system)
    : y(zeros(system.numElements())),
      z(zeros(system.numSpecies()))
    {
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

    // Return the stability indices of the phases
    auto phaseStabilityIndices(const EquilibriumState& state) const -> Vector
    {
        // The temperature and molar amounts of the species
        const double& T = state.temperature();
        const Vector& n = state.speciesAmounts();

        // Auxiliary variables
        const ChemicalSystem& system = state.system();
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

EquilibriumState::EquilibriumState()
: ChemicalState(), pimpl(new Impl())
{}

EquilibriumState::EquilibriumState(const ChemicalSystem& system)
: ChemicalState(system), pimpl(new Impl(system))
{}

EquilibriumState::EquilibriumState(const ChemicalState& state)
: ChemicalState(state), pimpl(new Impl(state.system()))
{}

EquilibriumState::EquilibriumState(const EquilibriumState& other)
: ChemicalState(other), pimpl(new Impl(*other.pimpl))
{}

EquilibriumState::~EquilibriumState()
{}

auto EquilibriumState::operator=(EquilibriumState other) -> EquilibriumState&
{
    ChemicalState::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumState::setSpeciesDualPotentials(const Vector& z) -> void
{
    pimpl->setSpeciesDualPotentials(z);
}

auto EquilibriumState::setElementDualPotentials(const Vector& y) -> void
{
    pimpl->setElementDualPotentials(y);
}

auto EquilibriumState::speciesDualPotentials() const -> const Vector&
{
    return pimpl->z;
}

auto EquilibriumState::elementDualPotentials() const -> const Vector&
{
    return pimpl->y;
}

auto EquilibriumState::phaseStabilityIndices() const -> Vector
{
    return pimpl->phaseStabilityIndices(*this);
}

auto EquilibriumState::output(std::string filename) -> void
{
    std::ofstream out(filename);
    out << *this;
}

auto operator<<(std::ostream& out, const EquilibriumState& state) -> std::ostream&
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
    const Vector molar_fractions = properties.molarFractions().val;
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
    const auto aqueous = properties.aqueous();
    const double I  = aqueous.ionicStrength().val;
    const double pH = aqueous.pH().val;
    const double pE = aqueous.pE().val;
    const double Eh = std::log(10)*R*T/F*pE;
    const double alk = aqueous.alkalinity().val;

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
    out << std::left << std::setw(25) << "Molar Fraction [mol/mol]";
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

} // namespace Reaktoro
