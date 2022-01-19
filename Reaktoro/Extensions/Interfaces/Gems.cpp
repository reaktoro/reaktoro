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

// #include "Gems.hpp"

// // C++ includes
// #include <map>
// #include <set>
// #include <vector>

// // Gems includes
// #define IPMGEMPLUGIN
// #define NOPARTICLEARRAY
// #include <gems/node.h>

// // Reaktoro includes
// #include <Reaktoro/Common/Constants.hpp>
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Common/TimeUtils.hpp>
// #include <Reaktoro/Core/ChemicalState.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>

// namespace Reaktoro {
// namespace {

// auto originalSpeciesName(const Gems& gems, unsigned index) -> std::string
// {
//     return gems.node()->pCSD()->DCNL[index];
// }

// auto uniqueSpeciesNames(const Gems& gems) -> std::vector<std::string>
// {
//     std::map<std::string, std::set<std::string>> species_names_in_phase;

//     unsigned offset = 0;
//     for(unsigned iphase = 0; iphase < gems.numPhases(); ++iphase)
//     {
//         const std::string phase_name = gems.phaseName(iphase);
//         for(unsigned i = 0; i < gems.numSpeciesInPhase(iphase); ++i)
//             species_names_in_phase[phase_name].insert(originalSpeciesName(gems, offset + i));
//         offset += gems.numSpeciesInPhase(iphase);
//     }

//     std::map<std::string, std::set<std::string>> species_found_in_phases;
//     for(unsigned i = 0; i < gems.numSpecies(); ++i)
//     {
//         std::string species_name = originalSpeciesName(gems, i);
//         for(const auto& pair : species_names_in_phase)
//         {
//             const auto& phase_name = pair.first;
//             const auto& species_names_in_this_phase = pair.second;
//             if(species_names_in_this_phase.count(species_name))
//                 species_found_in_phases[species_name].insert(phase_name);
//         }
//     }

//     std::vector<std::string> species_names;
//     species_names.reserve(gems.numSpecies());

//     for(unsigned i = 0; i < gems.numSpecies(); ++i)
//     {
//         std::string species_name = originalSpeciesName(gems, i);

//         if(species_found_in_phases[species_name].size() == 1)
//             species_names.push_back(species_name);
//         else
//         {
//             const Index iphase = gems.indexPhaseWithSpecies(i);
//             const std::string phase_name = gems.phaseName(iphase);
//             if(gems.numSpeciesInPhase(iphase) == 1)
//                 species_names.push_back(species_name);
//             else
//                 species_names.push_back(species_name + "(" + phase_name + ")");
//         }
//     }

//     return species_names;
// }

// } // namespace

// struct Gems::Impl
// {
//     /// The TNode instance from Gems
//     std::shared_ptr<TNode> node;

//     // The current temperature in GEMS TNode instance (in units of K)
//     double T;

//     // The current pressure in GEMS TNode instance (in units of Pa)
//     double P;

//     // The current molar amounts of all species in GEMS TNode instance (in units of mol)
//     VectorXr n;

//     /// The elapsed time of the equilibrate method (in units of s)
//     double elapsed_time = 0;

//     /// The options for Gems
//     GemsOptions options;

//     /// The unique names of the species
//     std::vector<std::string> species_names;

//     /// Construct a default Impl instance
//     Impl()
//     {}

//     /// Construct a default Impl instance
//     Impl(std::string filename)
//     {
//         // Initialize the GEMS `node` member
//         node = std::make_shared<TNode>();
//         if(node->GEM_init(filename.c_str()))
//             RuntimeError("Could not initialize the Gems object.",
//                 "Make sure the provided file exists relative to the working directory.");

//         //------------------------------------------------------------------------------------------------------
//         // The following parameters in GEMS have to be set to extremely small values to ensure that
//         // small molar amounts do not interfere with activity coefficient and chemical potential calculations
//         //------------------------------------------------------------------------------------------------------
//         // Reset the cutoff minimum amount of stable phase in GEMS (default: 1e-20)
//         node->pActiv()->GetActivityDataPtr()->DSM = 1e-300;

//         // Set the cutoff mole amount of water-solvent for aqueous phase elimination in GEMS (default: 1e-13)
//         node->pActiv()->GetActivityDataPtr()->XwMinM = 1e-300;

//         // Set the cutoff mole amount of solid sorbent for sorption phase elimination (default: 1e-13)
//         node->pActiv()->GetActivityDataPtr()->ScMinM = 1e-300;

//         // Set the cutoff mole amount for elimination of DC (species) in multi-component phase (default: 1e-33)
//         node->pActiv()->GetActivityDataPtr()->DcMinM = 1e-300;

//         // Set the cutoff mole amount for elimination of solution phases other than aqueous (default: 1e-20)
//         node->pActiv()->GetActivityDataPtr()->PhMinM = 1e-300;

//         // Set the cutoff effective molal ionic strength for calculation of aqueous activity coefficients (default: 1e-5)
//         node->pActiv()->GetActivityDataPtr()->ICmin = 1e-300;
//     }
// };

// Gems::Gems()
// : pimpl(new Impl())
// {}

// Gems::Gems(std::string filename)
// : pimpl(new Impl(filename))
// {
//     // Initialize the unique names of the species
//     pimpl->species_names = uniqueSpeciesNames(*this);

//     // Initialize the chemical state
//     set(temperature(), pressure(), speciesAmounts());
// }

// Gems::~Gems()
// {}

// auto Gems::temperature() const -> double
// {
//     return node()->Get_TK();
// }

// auto Gems::pressure() const -> double
// {
//     return node()->Get_P();
// }

// auto Gems::speciesAmounts() const -> VectorXr
// {
//     VectorXr n(numSpecies());
//     for(unsigned i = 0; i < n.size(); ++i)
//         n[i] = node()->Get_nDC(i);
//     return n;
// }

// auto Gems::numElements() const -> unsigned
// {
//     return node()->pCSD()->nIC;
// }

// auto Gems::numSpecies() const -> unsigned
// {
//     return node()->pCSD()->nDC;
// }

// auto Gems::numPhases() const -> unsigned
// {
//     return node()->pCSD()->nPH;
// }

// auto Gems::numSpeciesInPhase(Index iphase) const -> unsigned
// {
//     return node()->pCSD()->nDCinPH[iphase];
// }

// auto Gems::elementName(Index ielement) const -> std::string
// {
//     std::string name = node()->pCSD()->ICNL[ielement];
//     if(name == "Zz") name = "Z";
//     return name;
// }

// auto Gems::elementMolarMass(Index ielement) const -> double
// {
//     return node()->ICmm(ielement);
// }

// auto Gems::elementStoichiometry(Index ispecies, Index ielement) const -> double
// {
//     return node()->DCaJI(ispecies, ielement);
// }

// auto Gems::speciesName(Index ispecies) const -> std::string
// {
//     return pimpl->species_names[ispecies];
// }

// auto Gems::phaseName(Index iphase) const -> std::string
// {
//     return node()->pCSD()->PHNL[iphase];
// }

// auto Gems::properties(ThermoModelResult& res, double T, double P) -> void
// {
//     // Update the temperature and pressure of the Gems instance
//     set(T, P);

//     // The activity pointer from Gems
//     ACTIVITY* ap = node()->pActiv()->GetActivityDataPtr();

//     // The number of species and phases
//     const Index num_species = numSpecies();
//     const Index num_phases = numPhases();

//     // Set the thermodynamic properties of the species
//     for(Index i = 0; i < num_species; ++i)
//     {
//         res.standardPartialMolarGibbsEnergies().val[i] = node()->DC_G0(i, P, T, false);
//         res.standardPartialMolarEnthalpies().val[i] = node()->DC_H0(i, P, T);
//         res.standardPartialMolarVolumes().val[i] = node()->DC_V0(i, P, T);
//         res.standardPartialMolarHeatCapacitiesConstP().val[i] = node()->DC_Cp0(i, P, T);
//         res.standardPartialMolarHeatCapacitiesConstV().val[i] = node()->DC_Cp0(i, P, T);
//     }

//     Index offset = 0;
//     for(Index iphase = 0; iphase < num_phases; ++iphase)
//     {
//         // The number of species in the current phase
//         const Index size = numSpeciesInPhase(iphase);

//         // Set the ln activity constants of the species (non-zero for aqueous and gaseous species_
//         if(ap->PHC[iphase] == PH_AQUEL) // check if aqueous species
//         {
//             res.lnActivityConstants().val.segment(offset, size).fill(std::log(55.508472));
//             res.lnActivityConstants().val[ap->LO] = 0.0; // zero for water species
//         }
//         else if(ap->PHC[iphase] == PH_GASMIX) // check if gaseous species
//             res.lnActivityConstants().val.segment(offset, size).fill(std::log(1e-5 * P)); // ln(Pbar) for gases

//         offset += size;
//     }
// }

// auto Gems::properties(ChemicalModelResult& res, double T, double P, VectorXrConstRef n) -> void
// {
//     // Update the temperature, pressure, and species amounts of the Gems instance
//     set(T, P, n);

//     // The activity pointer from Gems
//     ACTIVITY* ap = node()->pActiv()->GetActivityDataPtr();

//     // The number of species and phases
//     const Index num_species = numSpecies();
//     const Index num_phases = numPhases();

//     // Set the ln activity coefficients and ln activities of the species in current phase
//     res.lnActivityCoefficients().val = VectorXr::Map(ap->lnGam, num_species);
//     res.lnActivities().val = VectorXr::Map(ap->lnAct, num_species);

//     // Set the molar derivatives of the activities
//     Index offset = 0;
//     for(Index iphase = 0; iphase < num_phases; ++iphase)
//     {
//         // The number of species in the current phase
//         const Index size = numSpeciesInPhase(iphase);

//         // The species amounts in the current phase
//         const auto np = n.segment(offset, size);

//         // Set the molar volume of current phase
//         res.phaseMolarVolumes().val[iphase] = (num_species == 1) ?
//             node()->DC_V0(offset, P, T) :
//             node()->Ph_Volume(iphase)/node()->Ph_Mole(iphase);

//         // Set d(ln(a))/dn to d(ln(x))/dn, where x is mole fractions
//         res.lnActivities().ddn.block(offset, offset, size, size) = -1.0/sum(np) * ones(size, size);
//         res.lnActivities().ddn.block(offset, offset, size, size).diagonal() += 1.0/np;

//         offset += size;
//     }
// }

// auto Gems::clone() const -> std::shared_ptr<Interface>
// {
//     return std::make_shared<Gems>(*this);
// }

// auto Gems::set(double T, double P) -> void
// {
//     pimpl->T = T;
//     pimpl->P = P;

//     node()->setTemperature(T);
//     node()->setPressure(P);
// }

// auto Gems::set(double T, double P, VectorXrConstRef n) -> void
// {
//     pimpl->T = T;
//     pimpl->P = P;
//     pimpl->n = n;

//     node()->setTemperature(T);
//     node()->setPressure(P);
//     node()->setSpeciation(n.data());

//     node()->updateStandardGibbsEnergies();
//     node()->initActivityCoefficients();
//     node()->updateConcentrations();
//     node()->updateActivityCoefficients();
//     node()->updateChemicalPotentials();
//     node()->updateActivities();
// }

// auto Gems::setOptions(const GemsOptions& options) -> void
// {
//     pimpl->options = options;
// }

// auto Gems::equilibrate(double T, double P, VectorXrConstRef b) -> void
// {
//     // Start timing
//     Time start = time();

//     // Set temperature and pressure
//     set(T, P);

//     // Set the molar amounts of the elements
//     for(unsigned i = 0; i < numElements(); ++i)
//         node()->pCNode()->bIC[i] = b[i];

//     // Solve the equilibrium problem with gems
//     node()->pCNode()->NodeStatusCH =
//         pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
//     node()->GEM_run(false);

//     // Finish timing
//     pimpl->elapsed_time = elapsed(start);
// }

// auto Gems::converged() const -> bool
// {
//     const auto status = node()->pCNode()->NodeStatusCH;
//     return status == OK_GEM_AIA || status == OK_GEM_SIA;
// }

// auto Gems::numIterations() const -> unsigned
// {
//     return node()->pCNode()->IterDone;
// }

// auto Gems::elapsedTime() const -> double
// {
//     return pimpl->elapsed_time;
// }

// auto Gems::node() const -> std::shared_ptr<TNode>
// {
//     return pimpl->node;
// }

// } // namespace Reaktoro
