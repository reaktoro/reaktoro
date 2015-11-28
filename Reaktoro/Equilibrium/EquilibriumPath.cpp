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

// C++ includes
#include <list>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalSensitivity.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Math/ODE.hpp>

namespace Reaktoro {

struct EquilibriumPath::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The options for the equilibrium path calculation and visualization
    EquilibriumPathOptions options;

    /// The partition of the chemical system
    Partition partition;

    /// The output instance of the equilibrium path calculation
    ChemicalOutput output;

    /// The plots of the equilibrium path calculation
    std::vector<ChemicalPlot> plots;

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
    auto solve(const ChemicalState& state_i, const ChemicalState& state_f) -> EquilibriumPathResult
    {
        // The result of this equilibrium path calculation
        EquilibriumPathResult result;

        // The number of equilibrium species
        const unsigned Ne = partition.numEquilibriumSpecies();

        // The indices of species and elements in the equilibrium partition
        const Indices& iequilibrium_species = partition.indicesEquilibriumSpecies();

        /// The temperatures at the initial and final chemical states
        const double T_i = state_i.temperature();
        const double T_f = state_f.temperature();

        /// The pressures at the initial and final chemical states
        const double P_i = state_i.pressure();
        const double P_f = state_f.pressure();

        /// The molar amounts of the elements in the equilibrium partition at the initial and final chemical states
        const Vector be_i = state_i.elementAmountsInSpecies(iequilibrium_species);
        const Vector be_f = state_f.elementAmountsInSpecies(iequilibrium_species);

        EquilibriumSolver equilibrium(system);
        equilibrium.setOptions(options.equilibrium);
        equilibrium.setPartition(partition);

        ChemicalState state = state_i;

        Vector dndT, dndP;
        Matrix dndb;

        ODEFunction f = [&](double t, const Vector& ne, Vector& res) -> int
        {
            const double T = T_i + t * (T_f - T_i);
            const double P = P_i + t * (P_f - P_i);
            const Vector be = be_i + t * (be_f - be_i);

            state.setTemperature(T);
            state.setPressure(P);

            result.equilibrium += equilibrium.solve(state, be);

            if(!result.equilibrium.optimum.succeeded) return 1;

            // The derivatives dn/dT, dn/dP, and dn/db for the equilibrium species
            dndT = rows(state.sensitivity().dndT, iequilibrium_species);
            dndP = rows(state.sensitivity().dndP, iequilibrium_species);
            dndb = rows(state.sensitivity().dndb, iequilibrium_species);

            // Calculate the right-hand side vector of the ODE
            res = dndT*(T_f - T_i) + dndP*(P_f - P_i) + dndb*(be_f - be_i);

            return 0;
        };

        ODEProblem problem;
        problem.setNumEquations(Ne);
        problem.setFunction(f);

        Vector ne = rows(state_i.speciesAmounts(), iequilibrium_species);

        // Adjust the absolute tolerance parameters for each component
        options.ode.abstols = options.ode.abstol * (ne + 1.0);

        // Ensure the iteration algorithm is not Newton
        options.ode.iteration = ODEIterationMode::Functional;

        ODESolver ode;
        ode.setOptions(options.ode);
        ode.setProblem(problem);

        double t = 0.0;

        ode.initialize(t, ne);

        // Initialize the output of the equilibrium path calculation
        if(output) output.open();

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        while(t < 1.0)
        {
            // Update the output with current state
            if(output) output.update(state, t);

            // Update the plots with current state
            for(auto& plot : plots) plot.update(state, t);

            // Integrate one time step only
            ode.integrate(t, ne);
        }

        // Update the output with the final state
        if(output) output.update(state_f, 1.0);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state_f, 1.0);

        return result;
    }
};

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

auto EquilibriumPath::setOptions(const EquilibriumPathOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto EquilibriumPath::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumPath::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumPath::solve(const ChemicalState& state_i, const ChemicalState& state_f) -> EquilibriumPathResult
{
    return pimpl->solve(state_i, state_f);
}

auto EquilibriumPath::output() -> ChemicalOutput
{
    pimpl->output = ChemicalOutput(pimpl->system);
    return pimpl->output;
}

auto EquilibriumPath::plot() -> ChemicalPlot
{
    pimpl->plots.push_back(ChemicalPlot(pimpl->system));
    return pimpl->plots.back();
}

auto EquilibriumPath::plots(unsigned num) -> std::vector<ChemicalPlot>
{
    for(unsigned i = 0; i < num; ++i) plot();
    return pimpl->plots;
}

} // namespace Reaktoro
