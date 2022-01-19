// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright © 2014-2022 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "EquilibriumReactions.hpp"

// // Reaktoro includes
// #include <Reaktoro/Common/Algorithms.hpp>
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>
// #include <Reaktoro/Core/Partition.hpp>
// #include <Reaktoro/Core/ReactionEquation.hpp>
// #include <Reaktoro/Core/Species.hpp>
// #include <Reaktoro/Core/ThermoProperties.hpp>
// #include <Reaktoro/Math/LU.hpp>
// #include <Reaktoro/Math/MathUtils.hpp>



// auto EquilibriumConditions::conservationMatrix() const -> MatrixXd
// {
//     const auto Nn = msystem.species().size();
//     const auto Ne = msystem.elements().size() + 1; // +1 = charge

//     const auto& inert_reactions = _details.restrictions.reactions_cannot_react;

//     const auto Nir = inert_reactions.size();

//     MatrixXd res = zeros(Ne + Nir, Nn);

//     auto fill_matrix_row = [=](const auto& pairs, auto row)
//     {
//         for(auto [ispecies, coeff] : pairs)
//             row[ispecies] = coeff;
//     };

//     auto i = 0;
//     for(const auto& pairs : inert_reactions)
//         fill_matrix_row(pairs, res.row(i++));

//     return res;
// }



// namespace Reaktoro {
// namespace {

// auto defaultMasterSpecies(const Partition& partition) -> Indices
// {
//     // The formula matrix of the equilibrium species
//     MatrixXdConstRef A = partition.formulaMatrixEquilibriumPartition();

//     // The number of elements and species in the equilibrium partition
//     const Index E = A.rows();
//     const Index N = A.cols();

//     // Return the number of species with a given element
//     auto num_species_with_element = [&](Index ielement)
//     {
//         Index count = 0;
//         for(Index j = 0; j < N; ++j)
//             if(A(ielement, j) != 0) ++count;
//         return count;
//     };

//     // Return the number of elements in a given species
//     auto num_elements_in_species = [&](Index ispecies)
//     {
//         Index count = 0;
//         for(Index i = 0; i < E; ++i)
//             if(A(i, ispecies) != 0) ++count;
//         return count;
//     };

//     // Return the number of species that have at least one of the elements in a given species
//     auto species_elemental_weight = [&](Index ispecies)
//     {
//         Index weight = 0;
//         for(Index i = 0; i < E; ++i)
//             if(A(i, ispecies) != 0) weight += num_species_with_element(i);
//         return weight;
//     };

//     // Return the indices of the species with a given element
//     auto indices_species_with_element = [&](Index ielement)
//     {
//         Indices indices;
//         for(Index j = 0; j < N; ++j)
//             if(A(ielement, j) != 0) indices.push_back(j);
//         return indices;
//     };

//     // Return for each element a potential master species
//     auto master_species_for_element = [&](Index ielement)
//     {
//         Indices ispecies = indices_species_with_element(ielement);
//         std::sort(ispecies.begin(), ispecies.end(),
//             [&](Index l, Index r) {
//                 if(num_elements_in_species(l) > num_elements_in_species(r))
//                     return false;
//                 if(num_elements_in_species(l) < num_elements_in_species(r))
//                     return true;
//                 if(species_elemental_weight(l) < species_elemental_weight(r))
//                     return false;
//                 if(species_elemental_weight(l) > species_elemental_weight(r))
//                     return true;
//                 if(sum(abs(A.col(l))) < sum(abs(A.col(r))))
//                     return true;
//                 return false;
//         });

//         return ispecies.front();
//     };

//     Indices imaster(E);
//     for(Index i = 0; i < E; ++i)
//         imaster[i] = master_species_for_element(i);

//     return imaster;
// }

// } // namespace

// struct EquilibriumReactions::Impl
// {
//     // The chemical system instance
//     ChemicalSystem system;

//     // The partition of the chemical species
//     Partition partition;

//     // The formula matrix of the equilibrium species
//     MatrixXd Ae;

//     // The indices of the equilibrium species
//     Indices iequilibrium;

//     // The number of species and elements in the equilibrium partition
//     Index Ne, Ee;

//     // The weighted formula matrix of the equilibrium species
//     MatrixXd We;

//     // The LU decomposition of the coefficient matrix `Abar`, where `P*Abar*Q = LU` and `Abar` is `Ae` without linearly dependent rows
//     LU lu;

//     // The indices of the master species
//     Indices imaster;

//     // The indices of the secondary species
//     Indices isecondary;

//     // The stoichiometric matrix of the equilibrium reactions
//     MatrixXd stoichiometric_matrix;

//     // The equations of the equilibrium reactions
//     std::vector<ReactionEquation> equations;

//     // Construct a Impl instance with given system
//     Impl(const ChemicalSystem& system)
//     : Impl(system, Partition(system))
//     {}

//     // Construct a Impl instance with given system and its partition
//     Impl(const ChemicalSystem& system, const Partition& partition)
//     : system(system), partition(partition)
//     {
//         // Initialize the formula matrix of the equilibrium species
//         Ae = partition.formulaMatrixEquilibriumPartition();

//         // Initialize the indices of the equilibrium species
//         iequilibrium = partition.indicesEquilibriumSpecies();

//         // The number of species and elements in the equilibrium partition
//         Ne = partition.numEquilibriumSpecies();
//         Ee = partition.numEquilibriumElements();

//         // Initialize the reactions with default master species
//         setMasterSpecies(defaultMasterSpecies(partition));
//     }

//     // Set the tentative master species with given priority of species to be selected as master species
//     auto setMasterSpecies(Indices imaster) -> void
//     {
//         // The number of equilibrium species
//         const Index num_species = iequilibrium.size();

//         // Initialize the weights of priority of the equilibrium species
//         VectorXr weights = ones(num_species);

//         // Give higher priority to the first indices in imaster
//         for(Index i = imaster.size(); i > 0; --i)
//             weights[imaster[i - 1]] = (imaster.size() - i + 1) * 100;

//         // Initialize the reactions with given weights
//         initialize(weights);
//     }

//     // Set the tentative master species with given names
//     auto setMasterSpecies(std::vector<std::string> species) -> void
//     {
//         // Get the global indices of the species
//         Indices ispecies = system.indicesSpecies(species);

//         // Convert the global indices to local indices (within the equilibrium partition)
//         for(Index& i : ispecies) i = index(iequilibrium, i);

//         // Assert all local indices are within bounds
//         for(Index i = 0; i < ispecies.size(); ++i)
//             Assert(ispecies[i] < iequilibrium.size(),
//                 "Could not initialize the equilibrium reactions with given master species.",
//                 "The master species `" + species[i] + "` is not present in the chemical system.");

//         // Finally set the master species
//         setMasterSpecies(ispecies);
//     }

//     // Initialize the equilibrium reactions based on weights of priority for master species
//     auto initialize(VectorXr weights) -> void
//     {
//         // Auxiliary references to the LU factors
//         const auto& U = lu.U;
//         const auto& Q = lu.Q;
//         const auto& r = lu.rank;

//         // Translate the weights so that the smallest value is one
//         weights.array() += 1 - min(weights);

//         // Compute the LU decomposition of Ae with scaling column-weights
//         lu.compute(Ae, weights);

//         // Initialize the indices of the master species
//         imaster = Indices(Q.indices().data(), Q.indices().data() + r);

//         // Initialize the indices of the secondary species
//         isecondary = Indices(Q.indices().data() + r, Q.indices().data() + Q.size());

//         // Convert local indices to global indices
//         for(Index& index : imaster) index = iequilibrium[index];
//         for(Index& index : isecondary) index = iequilibrium[index];

//         // Initialize the stoichiometric matrix
//         const Index num_primary = imaster.size();
//         const Index num_secondary = isecondary.size();
//         const Index num_species = num_primary + num_secondary;
//         const auto U1 = U.topLeftCorner(r, num_primary).triangularView<Eigen::Upper>();
//         const auto U2 = U.topRightCorner(r, num_secondary);
//         stoichiometric_matrix.resize(num_secondary, num_species);
//         stoichiometric_matrix.leftCols(num_primary) = tr(U1.solve(U2));
//         stoichiometric_matrix.rightCols(num_secondary) = -identity(num_secondary, num_secondary);
//         stoichiometric_matrix = stoichiometric_matrix * Q.inverse();

//         // Clean the stoichiometric matrix from round-off errors
//         cleanRationalNumbers(stoichiometric_matrix);

//         // Initialize the system of cannonical equilibrium reactions
//         equations.clear(); equations.reserve(num_secondary);
//         for(Index i = 0; i < num_secondary; ++i)
//         {
//             std::map<std::string, double> equation;
//             for(Index j = 0; j < num_species; ++j)
//                 if(stoichiometric_matrix(i, j) != 0)
//                     equation.emplace(system.species(iequilibrium[j]).name(), stoichiometric_matrix(i, j));
//             equations.push_back(ReactionEquation(equation));
//         }
//     }
// };

// EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system)
// : EquilibriumReactions(system, Partition(system))
// {
// }

// EquilibriumReactions::EquilibriumReactions(const ChemicalSystem& system, const Partition& partition)
// : pimpl(new Impl(system, partition))
// {

// }

// EquilibriumReactions::EquilibriumReactions(const EquilibriumReactions& other)
// : pimpl(new Impl(*other.pimpl))
// {
// }

// EquilibriumReactions::~EquilibriumReactions()
// {
// }

// auto EquilibriumReactions::operator=(EquilibriumReactions other) -> EquilibriumReactions&
// {
//     pimpl = std::move(other.pimpl);
//     return *this;
// }

// auto EquilibriumReactions::system() const -> const ChemicalSystem&
// {
//     return pimpl->system;
// }

// auto EquilibriumReactions::partition() const -> const Partition&
// {
//     return pimpl->partition;
// }


// auto EquilibriumReactions::setMasterSpecies(Indices ispecies) -> void
// {
//     pimpl->setMasterSpecies(ispecies);
// }

// auto EquilibriumReactions::setMasterSpecies(std::vector<std::string> species) -> void
// {
//     pimpl->setMasterSpecies(species);
// }

// auto EquilibriumReactions::indicesMasterSpecies() const -> Indices
// {
//     return pimpl->imaster;
// }

// auto EquilibriumReactions::indicesSecondarySpecies() const -> Indices
// {
//     return pimpl->isecondary;
// }

// auto EquilibriumReactions::equations() const -> std::vector<ReactionEquation>
// {
//     return pimpl->equations;
// }

// auto EquilibriumReactions::stoichiometricMatrix() const -> MatrixXd
// {
//     return pimpl->stoichiometric_matrix;
// }

// auto EquilibriumReactions::lu() const -> const LU&
// {
//     return pimpl->lu;
// }

// } // namespace Reaktoro
