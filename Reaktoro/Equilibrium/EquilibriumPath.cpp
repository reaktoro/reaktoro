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
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalQuantity.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
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

    /// The plots of the equilibrium path calculation
    std::vector<ChemicalPlot> plots;

    /// The chemical quantity instance
    ChemicalQuantity quantity;

    /// The collection of ChemicalState instances
    std::list<ChemicalState> states;

    /// The collection of `xi` values
    std::list<double> ts;

    /// Construct a default EquilibriumPath::Impl instance
    Impl()
    {}

    /// Construct a EquilibriumPath::Impl instance
    explicit Impl(const ChemicalSystem& system)
    : system(system), partition(system), quantity(system)
    {}

    /// Set the options for the equilibrium path calculation and visualization
    auto setOptions(const EquilibriumPathOptions& options_) -> void
    {
        options = options_;

        // Ensure the iteration algorithm is not Newton
        options.ode.iteration = ODEIterationMode::Functional;
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
        // The number of equilibrium species
        const unsigned Ne = partition.numEquilibriumSpecies();

        // The indices of equilibrium species
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

        Matrix dndb;

        ODEFunction f = [&](double xi, const Vector& ne, Vector& res) -> int
        {
            const double T = T_i + xi * (T_f - T_i);
            const double P = P_i + xi * (P_f - P_i);
            const Vector be = be_i + xi * (be_f - be_i);

            state.setTemperature(T);
            state.setPressure(P);

            auto result = equilibrium.solve(state, be);

            if(not result.optimum.succeeded) return 1;

//            todo Uncomment this once dn/dT and dn/dP can be calculated.
//            const Vector dndt = equilibrium.dndt(state);
//            const Vector dndp = equilibrium.dndp(state);
//            const Matrix dndb = equilibrium.dndb(state);
//            res = dndt*(T_f - T_i) + dndp*(P_f - P_i) + dndb*(be_f - be_i);

            dndb = equilibrium.dndb(state);

            res = dndb*(be_f - be_i);

            return 0;
        };

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        ODEProblem problem;
        problem.setNumEquations(Ne);
        problem.setFunction(f);

        Vector ne = rows(state_i.speciesAmounts(), iequilibrium_species);

        options.ode.abstols = options.ode.abstol * (ne + 1.0);

        ODESolver ode;
        ode.setOptions(options.ode);
        ode.setProblem(problem);

        double t = 0.0;

        ode.initialize(t, ne);

        while(t < 1.0)
        {
            // Update the plots with current state
            for(auto& plot : plots)
                plot.update(state, t);

            ts.push_back(t);
            states.push_back(state);
            ode.integrate(t, ne);
        }

        // Update the plots with final state
        for(auto& plot : plots)
            plot.update(state_f, 1.0);

        ts.push_back(1.0);
        states.push_back(state_f);
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

auto EquilibriumPath::solve(const ChemicalState& state_i, const ChemicalState& state_f) -> void
{
    pimpl->solve(state_i, state_f);
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
