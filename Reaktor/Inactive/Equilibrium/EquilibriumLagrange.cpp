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

#include "EquilibriumLagrange.hpp"

// C++ includes
#include <cmath>

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>
#include <Reaktor/Utils/SetUtils.hpp>

namespace Reaktor {

class EquilibriumLagrange::Impl
{
private:
    /// The definition of the chemical system
    ChemicalSystem system$;

    /// The partitioning of the chemical species in equilibrium, kinetic, and inert species
    Partitioning partitioning$;

    /// The Lagrange multipliers y of the equilibrium state
    Vector ye$;

    /// The Lagrange multipliers z of the equilibrium state
    Vector ze$;

    /// The stability indices of the equilibrium phases
    Vector stability_indices_equilibrium_phases$;

    /// The stability indices of all phases
    Vector stability_indices_phases$;

    /// The indices of the equilibrium species
    Indices idx_equilibrium_species$;

    /// The indices of the phases that contain equilibrium species
    Indices idx_equilibrium_phases$;

    /// The indices of the stable phases in the system (global indices)
    Indices idx_stable_phases$;

    /// The indices of the unstable phases in the system (global indices)
    Indices idx_unstable_phases$;

    /// The indices of the stable species in the system (global indices)
    Indices idx_stable_species$;

    /// The indices of the unstable species in the system (global indices)
    Indices idx_unstable_species$;

    /// The indices of the stable equilibrium phases in the system (local indices in the set of equilibrium species)
    Indices idx_equilibrium_stable_phases$;

    /// The indices of the unstable equilibrium phases in the system (local indices in the set of equilibrium species)
    Indices idx_equilibrium_unstable_phases$;

    /// The indices of the stable equilibrium species in the system (local indices in the set of equilibrium species)
    Indices idx_equilibrium_stable_species$;

    /// The indices of the unstable equilibrium species in the system (local indices in the set of equilibrium species)
    Indices idx_equilibrium_unstable_species$;

public:
    Impl(const ChemicalSystem& system)
    : Impl(system, Partitioning(system))
    {}

    Impl(const ChemicalSystem& system, const Partitioning& partitioning)
    : system$(system), partitioning$(partitioning),
      idx_equilibrium_species$(partitioning$.idxEquilibriumSpecies()),
      idx_equilibrium_phases$(partitioning$.idxPhasesWithEquilibriumSpecies())
    {}

    auto setMultipliers(const ChemicalState& state, const Vector& ye, const Vector& ze) -> void
    {
        ye$ = ye;
        ze$ = ze;
        initialiseStabilityIndices(state);
        initialiseStableUnstableIndices();
    }

    auto system() const -> const ChemicalSystem&
    {
        return system$;
    }

    auto partitioning() const -> const Partitioning&
    {
        return partitioning$;
    }

    auto lagrangeY() const -> const Vector&
    {
        return ye$;
    }

    auto lagrangeZ() const -> const Vector&
    {
        return ze$;
    }

    auto stabilityIndicesEquilibriumPhases() const -> const Vector&
    {
        return stability_indices_equilibrium_phases$;
    }

    auto stabilityIndicesPhases() const -> const Vector&
    {
        return stability_indices_phases$;
    }

    auto idxEquilibriumStablePhases() const -> Indices
    {
        return idx_equilibrium_stable_phases$;
    }

    auto idxEquilibriumUnstablePhases() const -> Indices
    {
        return idx_equilibrium_unstable_phases$;
    }

    auto idxEquilibriumStableSpecies() const -> Indices
    {
        return idx_equilibrium_stable_species$;
    }

    auto idxEquilibriumUnstableSpecies() const -> Indices
    {
        return idx_equilibrium_unstable_species$;
    }

    auto idxStablePhases() const -> Indices
    {
        return idx_stable_phases$;
    }

    auto idxUnstablePhases() const -> Indices
    {
        return idx_unstable_phases$;
    }

    auto idxStableSpecies() const -> Indices
    {
        return idx_stable_species$;
    }

    auto idxUnstableSpecies() const -> Indices
    {
        return idx_unstable_species$;
    }

    auto stablePhases() const -> std::vector<std::string>
    {
        return extract(system$.phasesNames(), idxStablePhases());
    }

    auto unstablePhases() const -> std::vector<std::string>
    {
        return extract(system$.phasesNames(), idxUnstablePhases());
    }

    auto stableSpecies() const -> std::vector<std::string>
    {
        return extract(system$.speciesNames(), idxStableSpecies());
    }

    auto unstableSpecies() const -> std::vector<std::string>
    {
        return extract(system$.speciesNames(), idxUnstableSpecies());
    }

    auto initialiseStabilityIndices(const ChemicalState& state) -> void
    {
        // The temperature times the universal gas constant
        const double RT = state.temperature() * universalGasConstant;

        // The molar fractions of all species
        const Vector x = state.molarFractions();

        // The extended Lagrange multipliers z to all species (not only equilibrium species)
        Vector z = zeros(system$.numSpecies());

        // Set the entries corresponding to the equilibrium species
        partitioning$.setEquilibriumRows(ze$, z);

        // Compute an auxiliary vector with the product wi' = xi * exp(-zi/RT)
        const Vector omega = x.array() / (z/RT).array().exp();

        // Initialise the stability indices of the equilibrium phases
        stability_indices_equilibrium_phases$ = zeros(system$.numPhases());

        // Iterate over all phases and compute the respective stability index
        for(unsigned i = 0; i < system$.numPhases(); ++i)
            stability_indices_equilibrium_phases$[i] = rows(system$.idxSpeciesInPhase(i), omega).sum();

        // Remove all stability indices of phases that do not have equilibrium species
        stability_indices_equilibrium_phases$ = rows(idx_equilibrium_phases$, stability_indices_equilibrium_phases$);

        // The natural log of 10
        const double ln10 = 2.30258509299;

        // Compute the last step of the stability indices of the equilibrium phases
        stability_indices_equilibrium_phases$ = stability_indices_equilibrium_phases$.array().log()/ln10;

        // Set the stability indices of all phases, with stability index zero for phase non-equilibrium phases
        stability_indices_phases$ = zeros(system$.numPhases());
        setRows(idx_equilibrium_phases$, stability_indices_equilibrium_phases$, stability_indices_phases$);
    }

    auto initialiseStableUnstableIndices() -> void
    {
        // Reset the indices of the stable and unstable species and phases
        idx_unstable_phases$.clear();
        idx_unstable_phases$.reserve(system$.numPhases());

        idx_stable_phases$.clear();
        idx_stable_phases$.reserve(system$.numPhases());

        idx_unstable_species$.clear();
        idx_unstable_species$.reserve(system$.numSpecies());

        idx_stable_species$.clear();
        idx_stable_species$.reserve(system$.numSpecies());

        idx_equilibrium_unstable_phases$.clear();
        idx_equilibrium_unstable_phases$.reserve(system$.numPhases());

        idx_equilibrium_stable_phases$.clear();
        idx_equilibrium_stable_phases$.reserve(system$.numPhases());

        idx_equilibrium_unstable_species$.clear();
        idx_equilibrium_unstable_species$.reserve(system$.numSpecies());

        idx_equilibrium_stable_species$.clear();
        idx_equilibrium_stable_species$.reserve(system$.numSpecies());

        // Collect the indices of the stable and unstable phases
        for(unsigned i = 0; i < idx_equilibrium_phases$.size(); ++i)
        {
            const Index idx_phase = idx_equilibrium_phases$[i];

            if(std::abs(stability_indices_equilibrium_phases$[i]) > 0.01)
            {
                idx_unstable_phases$.push_back(idx_phase);
                idx_equilibrium_unstable_phases$.push_back(i);
            }
            else
            {
                idx_stable_phases$.push_back(idx_phase);
                idx_equilibrium_stable_phases$.push_back(i);
            }
        }

        // Collect the indices of the stable and unstable species
        for(unsigned i = 0; i < idx_equilibrium_species$.size(); ++i)
        {
            const Index idx_species = idx_equilibrium_species$[i];
            const Index idx_phase = system$.idxPhaseWithSpecies(idx_species);

            if(contained(idx_phase, idx_stable_phases$))
            {
                idx_stable_species$.push_back(idx_species);
                idx_equilibrium_stable_species$.push_back(i);
            }
            else
            {
                idx_unstable_species$.push_back(idx_species);
                idx_equilibrium_unstable_species$.push_back(i);
            }
        }
    }
};

EquilibriumLagrange::EquilibriumLagrange(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumLagrange::EquilibriumLagrange(const ChemicalSystem& system, const Partitioning& partitioning)
: pimpl(new Impl(system, partitioning))
{}

EquilibriumLagrange::EquilibriumLagrange(const EquilibriumLagrange& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumLagrange::~EquilibriumLagrange()
{}

auto EquilibriumLagrange::operator=(EquilibriumLagrange other) -> EquilibriumLagrange&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumLagrange::setMultipliers(const ChemicalState& state, const Vector& ye, const Vector& ze) -> void
{
    return pimpl->setMultipliers(state, ye, ze);
}

auto EquilibriumLagrange::system() const -> const ChemicalSystem&
{
    return pimpl->system();
}

auto EquilibriumLagrange::partitioning() const -> const Partitioning&
{
    return pimpl->partitioning();
}

auto EquilibriumLagrange::lagrangeY() const -> const Vector&
{
    return pimpl->lagrangeY();
}

auto EquilibriumLagrange::lagrangeZ() const -> const Vector&
{
    return pimpl->lagrangeZ();
}

auto EquilibriumLagrange::stabilityIndicesEquilibriumPhases() const -> const Vector&
{
    return pimpl->stabilityIndicesEquilibriumPhases();
}

auto EquilibriumLagrange::stabilityIndicesPhases() const -> const Vector&
{
    return pimpl->stabilityIndicesPhases();
}

auto EquilibriumLagrange::idxEquilibriumStablePhases() const -> Indices
{
    return pimpl->idxEquilibriumStablePhases();
}

auto EquilibriumLagrange::idxEquilibriumUnstablePhases() const -> Indices
{
    return pimpl->idxEquilibriumUnstablePhases();
}

auto EquilibriumLagrange::idxEquilibriumStableSpecies() const -> Indices
{
    return pimpl->idxEquilibriumStableSpecies();
}

auto EquilibriumLagrange::idxEquilibriumUnstableSpecies() const -> Indices
{
    return pimpl->idxEquilibriumUnstableSpecies();
}

auto EquilibriumLagrange::idxStablePhases() const -> Indices
{
    return pimpl->idxStablePhases();
}

auto EquilibriumLagrange::idxUnstablePhases() const -> Indices
{
    return pimpl->idxUnstablePhases();
}

auto EquilibriumLagrange::idxStableSpecies() const -> Indices
{
    return pimpl->idxStableSpecies();
}

auto EquilibriumLagrange::idxUnstableSpecies() const -> Indices
{
    return pimpl->idxUnstableSpecies();
}

auto EquilibriumLagrange::stablePhases() const -> std::vector<std::string>
{
    return pimpl->stablePhases();
}

auto EquilibriumLagrange::unstablePhases() const -> std::vector<std::string>
{
    return pimpl->unstablePhases();
}

auto EquilibriumLagrange::stableSpecies() const -> std::vector<std::string>
{
    return pimpl->stableSpecies();
}

auto EquilibriumLagrange::unstableSpecies() const -> std::vector<std::string>
{
    return pimpl->unstableSpecies();
}

} /* namespace Reaktor */
