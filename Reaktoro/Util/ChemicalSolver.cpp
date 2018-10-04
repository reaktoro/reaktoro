// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.
//
//#include "ChemicalSolver.hpp"
//
//// C++ includes
//#include <vector>
//
//// Reaktoro includes
//#include <Reaktoro/Common/Exception.hpp>
//#include <Reaktoro/Core/ChemicalProperties.hpp>
//#include <Reaktoro/Core/ChemicalState.hpp>
//#include <Reaktoro/Core/ChemicalSystem.hpp>
//#include <Reaktoro/Core/Partition.hpp>
//#include <Reaktoro/Core/ReactionSystem.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
//#include <Reaktoro/Kinetics/KineticSolver.hpp>
//#include <Reaktoro/Util/ChemicalField.hpp>
//
//namespace Reaktoro {
//
//struct ChemicalSolver::Impl
//{
//    /// The chemical system instance
//    ChemicalSystem system;
//
//    /// The reaction system instance
//    ReactionSystem reactions;
//
//    /// The number of field points
//    Index npoints;
//
//    /// The partitioning of the chemical system
//    Partition partition;
//
//    /// The number of species and elements in the system
//    Index N, E;
//
//    /// The number of species and elements in the equilibrium partition
//    Index Ne, Ee, Nk, Nfp;
//
//    /// The number of components
//    Index Nc;
//
//    /// The chemical states at each point in the field
//    std::vector<ChemicalState> states;
//
//    /// The chemical properties at each point in the field
//    std::vector<ChemicalProperties> properties;
//
//    /// The equilibrium solver
//    EquilibriumSolver equilibriumsolver;
//
//    /// The kinetic solver
//    KineticSolver kineticsolver;
//
//    /// The equilibrium sensitivity at every field point
//    std::vector<EquilibriumSensitivity> sensitivities;
//
//    /// The molar amounts of the chemical components at every field point
//    std::vector<Vector> c;
//
//    /// The kinetic rates of the chemical components and their derivatives at every field point (in units of mol/s)
//    std::vector<ChemicalField> rc;
//
//    /// The molar amounts of equilibrium species and their derivatives at every field point
//    std::vector<ChemicalField> ne;
//
//    /// The porosity at every field point and their derivatives
//    ChemicalField porosity;
//
//    /// The saturations of the fluid phases and their derivatives at every field point
//    std::vector<ChemicalField> fluid_saturations;
//
//    /// The densities of the fluid phases and their derivatives at every field point (in units of kg/m3)
//    std::vector<ChemicalField> fluid_densities;
//
//    /// The volumes of the fluid phases and their derivatives at every field point (in units of m3)
//    std::vector<ChemicalField> fluid_volumes;
//
//    /// The total volume of the fluid phases and their derivatives at every field point (in units of m3)
//    ChemicalField fluid_total_volume;
//
//    /// The total volume of the solid phases and their derivatives at every field point (in units of m3)
//    ChemicalField solid_total_volume;
//
//    /// Construct a default Impl instance
//    Impl()
//    {}
//
//    /// Construct a custom Impl instance with given chemical system
//    Impl(const ChemicalSystem& system, Index npoints)
//    : system(system),
//      npoints(npoints),
//      states(npoints, ChemicalState(system)),
//      properties(npoints),
//      equilibriumsolver(system)
//    {
//        setPartition(Partition(system));
//    }
//
//    /// Construct a custom Impl instance with given reaction system
//    Impl(const ReactionSystem& reactions, Index npoints)
//    : system(reactions.system()),
//      reactions(reactions),
//      npoints(npoints),
//      states(npoints, ChemicalState(system)),
//      properties(npoints),
//      equilibriumsolver(system),
//      kineticsolver(reactions)
//    {
//        // Initialize the number of species and elements in the system
//        N = system.numSpecies();
//        E = system.numElements();
//
//        // Initialize the default partition of the chemical system
//        setPartition(Partition(system));
//    }
//
//    /// Set the partition of the chemical system
//    auto setPartition(const Partition& partition_) -> void
//    {
//        // Set the partition of the chemical solver
//        partition = partition_;
//
//        // Initialize the number-type variables
//        Ne  = partition.numEquilibriumSpecies();
//        Nk  = partition.numKineticSpecies();
//        Nfp = partition.numFluidPhases();
//        Ee  = partition.numEquilibriumElements();
//        Nc  = Ee + Nk;
//
//        // Set the partition of the equilibrium and kinetic solvers
//        if(Ne) equilibriumsolver.setPartition(partition);
//        if(Nk) kineticsolver.setPartition(partition);
//
//        // Initialize the sensitivities member
//        sensitivities.resize(npoints);
//    }
//
//    /// Equilibrate the chemical state at every field point.
//    auto equilibrate(Array<double> T, Array<double> P, Array<double> b) -> void
//    {
//        Assert(T.size == npoints,
//            "Could not perform equilibrium calculations.",
//            "Expecting the same number of temperature values as there are field points.");
//
//        Assert(P.size == npoints,
//            "Could not perform equilibrium calculations.",
//            "Expecting the same number of pressure values as there are field points.");
//
//        Assert(b.size == npoints * Ee,
//            "Could not perform equilibrium calculations.",
//            "Expecting, for each equilibrium element, the same number of amount "
//            "values as there are field points.");
//
//        for(Index k = 0; k < npoints; ++k)
//        {
//            const auto Tk = T.data[k];
//            const auto Pk = P.data[k];
//            const auto bk = b.data + k*Ee;
//            equilibriumsolver.solve(states[k], Tk, Pk, bk);
//            properties[k] = states[k].properties();
//            sensitivities[k] = equilibriumsolver.sensitivity();
//        }
//    }
//
//    /// Equilibrate the chemical state at every field point.
//    auto equilibrate(Array<double> T, Array<double> P, Grid<double> b) -> void
//    {
//        Assert(T.size == npoints,
//            "Could not perform equilibrium calculations.",
//            "Expecting the same number of temperature values as there are field points.");
//
//        Assert(P.size == npoints,
//            "Could not perform equilibrium calculations.",
//            "Expecting the same number of pressure values as there are field points.");
//
//        Assert(b.rows == Ee && b.cols == npoints,
//            "Could not perform equilibrium calculations.",
//            "Expecting, for each equilibrium element, the same number of amount "
//            "values as there are field points.");
//
//        Vector bk(Ee);
//        for(Index k = 0; k < npoints; ++k)
//        {
//            const auto Tk  = T.data[k];
//            const auto Pk  = P.data[k];
//            for(Index j = 0; j < Ee; ++j)
//                bk[j] = b.data[j][k];
//            equilibriumsolver.solve(states[k], Tk, Pk, bk);
//            properties[k] = states[k].properties();
//            sensitivities[k] = equilibriumsolver.sensitivity();
//        }
//    }
//
//    /// React the chemical state at every field point.
//    auto react(double t, double dt) -> void
//    {
//        for(Index k = 0; k < npoints; ++k)
//        {
//            kineticsolver.solve(states[k], t, dt);
//            properties[k] = states[k].properties();
//        }
//    }
//
//    /// Update the molar amounts of the chemical components at every field point.
//    auto updateComponentAmounts() -> void
//    {
//        // Check if member `c` is initialized
//        if(c.size() != Nc)
//            c.resize(Nc, Vector(npoints));
//
//        // The indices of the equilibrium elements and kinetic species
//        const Indices& ies = partition.indicesEquilibriumSpecies();
//        const Indices& iee = partition.indicesEquilibriumElements();
//        const Indices& iks = partition.indicesKineticSpecies();
//
//        // The molar amounts of elements in the equilibrium species
//        Vector bes(E);
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Get the molar amounts of all species
//            VectorConstRef n = states[k].speciesAmounts();
//
//            // Calculate the molar amounts of the elements in the equilibrium species
//            bes = states[k].elementAmountsInSpecies(ies);
//
//            // Loop over all equilibrium elements
//            for(Index j = 0; j < Ee; ++j)
//                c[j][k] = bes[iee[j]];
//
//            // Loop over all kinetic species
//            for(Index i = 0; i < Nk; ++i)
//                c[i + Ee][k] = n[iks[i]];
//        }
//    }
//
//    /// Update the molar amounts of the equilibrium species and their derivatives at every field point.
//    auto updateEquilibriumSpeciesAmounts() -> void
//    {
//        // Check if member `ne` is initialized
//        if(ne.size() != Ne)
//            ne.resize(Ne, ChemicalField(partition, npoints));
//
//        // The indices of the equilibrium species
//        const Indices& ies = partition.indicesEquilibriumSpecies();
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Loop over all equilibrium species
//            for(Index i = 0; i < Ne; ++i)
//            {
//                // Set the molar amount of the current equilibrium species
//                ne[i].val()[k] = states[k].speciesAmount(ies[i]);
//
//                // Set the sensitivity of the current equilibrium species w.r.t. temperature
//                ne[i].ddT()[k] = sensitivities[k].dnedT[i];
//
//                // Set the sensitivity of the current equilibrium species w.r.t. pressure
//                ne[i].ddP()[k] = sensitivities[k].dnedP[i];
//
//                // Set the sensitivity of the current equilibrium species w.r.t. molar amount of each equilibrium element
//                for(Index j = 0; j < Ee; ++j)
//                    ne[i].ddbe()[j][k] = sensitivities[k].dnedbe(i, j);
//            }
//        }
//    }
//
//    /// Update the porosity and their derivatives at every field point.
//    auto updatePorosity() -> void
//    {
//        // Check if member `porosity` is initialized
//        if(!porosity.size())
//            porosity = ChemicalField(partition, npoints);
//
//        // The porosity
//        ChemicalScalar phi;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Calculate porosity at current field point
//            phi = 1.0 - properties[k].solidVolume();
//            porosity.set(k, phi, sensitivities[k]);
//        }
//    }
//
//    /// Get the saturation of each fluid phase at every field point.
//    auto updateFluidSaturations() -> void
//    {
//        // Check if member `fluid_saturations` is initialized
//        if(fluid_saturations.size() != Nfp)
//            fluid_saturations.resize(Nfp, ChemicalField(partition, npoints));
//
//        // The indices of the fluid phases
//        const Indices& ifp = partition.indicesFluidPhases();
//
//        // The volumes of the fluid phases and their saturations
//        ChemicalVector fluid_volumes, fluid_saturation;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Get the volumes of the fluid phases
//            fluid_volumes = rows(properties[k].phaseVolumes(), ifp);
//
//            // Compute the saturation of all fluid phases
//            fluid_saturation = fluid_volumes/sum(fluid_volumes);
//
//            // Loop over all fluid phases
//            for(Index j = 0; j < Nfp; ++j)
//                fluid_saturations[j].set(k, fluid_saturation[j], sensitivities[k]);
//        }
//    }
//
//    /// Get the density of each fluid phase at every field point.
//    auto updateFluidDensities() -> void
//    {
//        // Check if member `fluid_densities` is initialized
//        if(fluid_densities.size() != Nfp)
//            fluid_densities.resize(Nfp, ChemicalField(partition, npoints));
//
//        // The indices of the fluid phases
//        const Indices& ifp = partition.indicesFluidPhases();
//
//        // The densities of all phases at a field point
//        ChemicalVector rho;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Get the densities of all phases
//            rho = rows(properties[k].phaseDensities(), ifp);
//
//            // Loop over all fluid phases
//            for(Index j = 0; j < Nfp; ++j)
//                fluid_densities[j].set(k, rho[j], sensitivities[k]);
//        }
//    }
//
//    /// Update the volumes of each fluid phase at every field point.
//    auto updateFluidVolumes() -> void
//    {
//        // Check if member `fluid_volumes` is initialized
//        if(fluid_volumes.size() != Nfp)
//            fluid_volumes.resize(Nfp, ChemicalField(partition, npoints));
//
//        // The indices of the fluid phases
//        const Indices& ifp = partition.indicesFluidPhases();
//
//        // The volumes of all fluid phases
//        ChemicalVector volumes;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Get the volumes of all phases
//            volumes = rows(properties[k].phaseVolumes(), ifp);
//
//            // Loop over all fluid phases
//            for(Index j = 0; j < Nfp; ++j)
//                fluid_volumes[j].set(k, volumes[j], sensitivities[k]);
//        }
//    }
//
//    /// Update the total volume of the fluid phases at every field point.
//    auto updateFluidTotalVolume() -> void
//    {
//        // Check if member `fluid_total_volume` is initialized
//        if(!fluid_total_volume.size())
//            fluid_total_volume = ChemicalField(partition, npoints);
//
//        // The indices of the fluid phases
//        const Indices& ifp = partition.indicesFluidPhases();
//
//        // The total fluid volume at a field point
//        ChemicalScalar total_volume;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Calculate the total fluid volume at current field point
//            total_volume = sum(rows(properties[k].phaseVolumes(), ifp));
//            fluid_total_volume.set(k, total_volume, sensitivities[k]);
//        }
//    }
//
//    /// Update the total volume of the solid phases at every field point.
//    auto updateSolidTotalVolume() -> void
//    {
//        // Check if member `solid_total_volume` is initialized
//        if(!solid_total_volume.size())
//            solid_total_volume = ChemicalField(partition, npoints);
//
//        // The indices of the solid phases
//        const Indices& isp = partition.indicesSolidPhases();
//
//        // The total solid volume at a field point
//        ChemicalScalar total_volume;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            // Calculate the total solid volume at current field point
//            total_volume = sum(rows(properties[k].phaseVolumes(), isp));
//            solid_total_volume.set(k, total_volume, sensitivities[k]);
//        }
//    }
//
//    /// Update the kinetic rates of the chemical components at every field point.
//    auto updateComponentRates() -> void
//    {
//        // Check if member `rc` is initialized
//        if(rc.size() != Nc)
//            rc.resize(Nc, ChemicalField(partition, npoints));
//
//        // Skip the calculation if there are no reactions
//        if(reactions.numReactions() == 0)
//            return;
//
//        // The indices of the elements and species in the equilibrium partition
//        const Indices& iee = partition.indicesEquilibriumElements();
//        const Indices& ies = partition.indicesEquilibriumSpecies();
//
//        // The indices of the species in the kinetic partition
//        const Indices& iks = partition.indicesKineticSpecies();
//
//        // The formula matrix w.r.t. the species and elements in the equilibrium partition
//        const Matrix We = submatrix(system.formulaMatrix(), iee, ies);
//
//        // The stoichiometric matrices w.r.t. the equilibrium and kinetic species
//        const Matrix Se = cols(reactions.stoichiometricMatrix(), ies);
//        const Matrix Sk = cols(reactions.stoichiometricMatrix(), iks);
//
//        // Initialise the coefficient matrix `A` of the chemical kinetics problem
//        Matrix A(Ee + Nk, reactions.numReactions());
//        A.topRows(Ee) = We * tr(Se);
//        A.bottomRows(Nk) = tr(Sk);
//
//        ChemicalVector rates, r;
//
//        // Loop over all field points
//        for(Index k = 0; k < npoints; ++k)
//        {
//            r = reactions.rates(properties[k]);
//
//            rates.val = A * r.val;
//            rates.ddT = A * r.ddT;
//            rates.ddP = A * r.ddP;
//            rates.ddn = A * r.ddn;
//
//            // Loop over all equilibrium elements
//            for(Index j = 0; j < Ee; ++j)
//                rc[j].set(k, rates[j], sensitivities[k]);
//
//            // Loop over all kinetic species
//            for(Index i = 0; i < Nk; ++i)
//                rc[i + Ee].set(k, rates[i + Ee], sensitivities[k]);
//        }
//    }
//};
//
//ChemicalSolver::ChemicalSolver()
//: pimpl(new Impl())
//{}
//
//ChemicalSolver::ChemicalSolver(const ChemicalSystem& system, Index npoints)
//: pimpl(new Impl(system, npoints))
//{}
//
//ChemicalSolver::ChemicalSolver(const ReactionSystem& reactions, Index npoints)
//: pimpl(new Impl(reactions, npoints))
//{}
//
//auto ChemicalSolver::numPoints() const -> Index
//{
//    return pimpl->npoints;
//}
//
//auto ChemicalSolver::numEquilibriumElements() const -> Index
//{
//    return pimpl->Ee;
//}
//
//auto ChemicalSolver::numKineticSpecies() const -> Index
//{
//    return pimpl->Nk;
//}
//
//auto ChemicalSolver::numComponents() const -> Index
//{
//    return pimpl->Nc;
//}
//
//auto ChemicalSolver::setPartition(const Partition& partition) -> void
//{
//    pimpl->setPartition(partition);
//}
//
//auto ChemicalSolver::setStates(const ChemicalState& state) -> void
//{
//    for(Index k = 0; k < pimpl->npoints; ++k)
//        pimpl->states[k] = state;
//}
//
//auto ChemicalSolver::setStates(const Array<ChemicalState>& states) -> void
//{
//    Assert(states.size == pimpl->npoints,
//        "Could not set the chemical states at every field point.",
//        "Expecting the same number of chemical states as there are field points.");
//    for(Index k = 0; k < states.size; ++k)
//        pimpl->states[k] = states.data[k];
//}
//
//auto ChemicalSolver::setStateAt(Index ipoint, const ChemicalState& state) -> void
//{
//    Assert(ipoint < pimpl->npoints,
//        "Could not set the chemical state at given field point.",
//        "Expecting a field point index smaller than the number of field points.");
//    pimpl->states[ipoint] = state;
//}
//
//auto ChemicalSolver::setStateAt(const Array<Index>& ipoints, const ChemicalState& state) -> void
//{
//    Assert(ipoints.size < pimpl->npoints,
//        "Could not set the chemical state at given field points.",
//        "Expecting field point indices smaller than the number of field points.");
//    for(Index k = 0; k < ipoints.size; ++k)
//        setStateAt(ipoints.data[k], state);
//}
//
//auto ChemicalSolver::setStateAt(const Array<Index>& ipoints, const Array<ChemicalState>& states) -> void
//{
//    Assert(ipoints.size < pimpl->npoints,
//        "Could not set the chemical state at given field points.",
//        "Expecting field point indices smaller than the number of field points.");
//    Assert(ipoints.size == states.size,
//        "Could not set the chemical state at given field points.",
//        "Expecting the same number of field point indices and chemical states.");
//    for(Index k = 0; k < ipoints.size; ++k)
//        setStateAt(ipoints.data[k], states.data[k]);
//}
//
//auto ChemicalSolver::equilibrate(Array<double> T, Array<double> P, Array<double> be) -> void
//{
//    pimpl->equilibrate(T, P, be);
//}
//
//auto ChemicalSolver::equilibrate(Array<double> T, Array<double> P, Grid<double> be) -> void
//{
//    pimpl->equilibrate(T, P, be);
//}
//
//auto ChemicalSolver::react(double t, double dt) -> void
//{
//    pimpl->react(t, dt);
//}
//
//auto ChemicalSolver::state(Index i) const -> const ChemicalState&
//{
//    return pimpl->states[i];
//}
//
//auto ChemicalSolver::states() const -> const std::vector<ChemicalState>&
//{
//    return pimpl->states;
//}
//
//auto ChemicalSolver::componentAmounts() -> const std::vector<Vector>&
//{
//    pimpl->updateComponentAmounts();
//    return pimpl->c;
//}
//
//auto ChemicalSolver::equilibriumSpeciesAmounts() -> const std::vector<ChemicalField>&
//{
//    pimpl->updateEquilibriumSpeciesAmounts();
//    return pimpl->ne;
//}
//
//auto ChemicalSolver::porosity() -> const ChemicalField&
//{
//    pimpl->updatePorosity();
//    return pimpl->porosity;
//}
//
//auto ChemicalSolver::fluidSaturations() -> const std::vector<ChemicalField>&
//{
//    pimpl->updateFluidSaturations();
//    return pimpl->fluid_saturations;
//}
//
//auto ChemicalSolver::fluidDensities() -> const std::vector<ChemicalField>&
//{
//    pimpl->updateFluidDensities();
//    return pimpl->fluid_densities;
//}
//
//auto ChemicalSolver::fluidVolumes() -> const std::vector<ChemicalField>&
//{
//    pimpl->updateFluidVolumes();
//    return pimpl->fluid_volumes;
//}
//
//auto ChemicalSolver::fluidTotalVolume() -> const ChemicalField&
//{
//    pimpl->updateFluidTotalVolume();
//    return pimpl->fluid_total_volume;
//}
//
//auto ChemicalSolver::solidTotalVolume() -> const ChemicalField&
//{
//    pimpl->updateSolidTotalVolume();
//    return pimpl->solid_total_volume;
//}
//
//auto ChemicalSolver::componentRates() -> const std::vector<ChemicalField>&
//{
//    pimpl->updateComponentRates();
//    return pimpl->rc;
//}
//
//}  // namespace Reaktoro
