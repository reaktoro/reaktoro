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

#include "EquilibriumPath.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {

struct EquilibriumPath::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The options for the equilibrium path calculation and visualization
    EquilibriumPathOptions options;

    /// The partition of the chemical system
    Partition partition;

    /// The vector of ChemicalState instances
    std::vector<ChemicalState> states;

    /// Construct a default EquilibriumPath::Impl instance
    Impl()
    {}

    /// Construct a EquilibriumPath::Impl instance
    explicit Impl(const ChemicalSystem& system)
    : system(system), partition(system)
    {}

    /// Set the options for the equilibrium path calculation and visualization
    auto setOptions(const EquilibriumPathOptions& options_) -> void
    {
        options = options_;
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;
    }

    /// Set the partition of the chemical system using a formatted string
    auto setPartition(std::string partition) -> void
    {
        setPartition(Partition(system, partition));
    }

    /// Solve the path of equilibrium states between two chemical states
    auto solve(const ChemicalState& state_i, const ChemicalState& state_f) -> void
    {
        // The indices of equilibrium species
        auto iequilibrium_species = partition.indicesEquilibriumSpecies();

        // Resize the collection of chemical states
        states.resize(options.num_points);

        // Set the first and last chemical state
        states.front() = state_i;
        states.back() = state_f;

        /// The temperatures at the initial and final chemical states
        const double T_i = state_i.temperature();
        const double T_f = state_f.temperature();

        /// The pressures at the initial and final chemical states
        const double P_i = state_i.pressure();
        const double P_f = state_f.pressure();

        /// The molar amounts of the elements in the equilibrium partition at the initial and final chemical states
        const Vector be_i = state_i.elementAmountsInSpecies(iequilibrium_species);
        const Vector be_f = state_f.elementAmountsInSpecies(iequilibrium_species);

        EquilibriumSolver solver(system);
        solver.setOptions(options.equilibrium);
        solver.setPartition(partition);

        for(unsigned i = 1; i < states.size() - 1; ++i)
        {
            const double alpha = static_cast<double>(i)/(states.size() - 1);

            const double T = T_i + alpha * (T_f - T_i);
            const double P = P_i + alpha * (P_f - P_i);

            states[i] = states[i - 1];
            states[i].setTemperature(T);
            states[i].setPressure(P);

            const Vector be = be_i + alpha * (be_f - be_i);

            solver.solve(states[i], be);
        }
    }

    /// Output the equilibrium path
    auto output(std::string list) -> void
    {
        auto quantities = split(list, " ");

        // Output the header
        for(std::string quantity : quantities)
            std::cout << std::left << std::setw(20) << quantity;
        std::cout << std::endl;

        // Output the quantities for each chemical state
        for(const ChemicalState& state : states)
        {
            for(std::string quantity : quantities)
                std::cout << std::left << std::setw(20) << extract(state, quantity);
            std::cout << std::endl;
        }
    }

    /// Plot the equilibrium path
    auto plot(std::string list) -> void
    {

    }
};

EquilibriumPath::EquilibriumPath()
: pimpl(new Impl())
{}

EquilibriumPath::EquilibriumPath(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumPath::EquilibriumPath(const EquilibriumPath& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumPath::~EquilibriumPath()
{}

auto EquilibriumPath::operator=(EquilibriumPath other) -> EquilibriumPath&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumPath::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumPath::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumPath::solve(const ChemicalState& state_i, const ChemicalState& state_f) -> void
{
    pimpl->solve(state_i, state_f);
}

auto EquilibriumPath::output(std::string list) -> void
{
    pimpl->output(list);
}

auto EquilibriumPath::plot(std::string list) -> void
{
    pimpl->plot(list);
}

} // namespace Reaktoro
