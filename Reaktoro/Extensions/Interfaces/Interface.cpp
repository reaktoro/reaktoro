// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2021 Allan Leal
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

// #include "Interface.hpp"

// // C++ includes
// #include <map>
// #include <vector>

// // Reaktoro includes
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Core/ChemicalState.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>
// #include <Reaktoro/Core/Element.hpp>
// #include <Reaktoro/Core/Phase.hpp>
// #include <Reaktoro/Core/Species.hpp>

// namespace Reaktoro {
// namespace {

// /// Return the Element instances and their stoichiometries that compose a species
// auto elementsInSpecies(const Interface& interface, std::vector<Element> elements, Index ispecies) -> std::map<Element, double>
// {
//     std::map<Element, double> res;
//     for(unsigned i = 0; i < interface.numElements(); ++i)
//         if(interface.elementStoichiometry(ispecies, i) != 0.0)
//             res.emplace(elements[i], interface.elementStoichiometry(ispecies, i));
//     return res;
// }

// /// Return the Species instances that compose a phase
// auto speciesInPhase(const Interface& interface, std::vector<Species> species, Index iphase) -> std::vector<Species>
// {
//     const unsigned ifirst = interface.indexFirstSpeciesInPhase(iphase);
//     const unsigned nspecies = interface.numSpeciesInPhase(iphase);
//     const auto begin = species.begin() + ifirst;
//     const auto end = begin + nspecies;
//     return std::vector<Species>(begin, end);
// }

// } // namespace

// Interface::~Interface()
// {}

// auto Interface::formulaMatrix() const -> MatrixXd
// {
//     const unsigned E = numElements();
//     const unsigned N = numSpecies();
//     MatrixXd A(E, N);
//     for(unsigned i = 0; i < N; ++i)
//         for(unsigned j = 0; j < E; ++j)
//             A(j, i) = elementStoichiometry(i, j);
//     return A;
// }

// auto Interface::indexElement(std::string element) const -> Index
// {
//     const Index size = numElements();
//     for(unsigned i = 0; i < size; ++i)
//         if(elementName(i) == element)
//             return i;
//     return size;
// }

// auto Interface::indexSpecies(std::string species) const -> Index
// {
//     const Index size = numSpecies();
//     for(unsigned i = 0; i < size; ++i)
//         if(speciesName(i) == species)
//             return i;
//     return size;
// }

// auto Interface::indexPhase(std::string phase) const -> Index
// {
//     const Index size = numPhases();
//     for(unsigned i = 0; i < size; ++i)
//         if(phaseName(i) == phase)
//             return i;
//     return size;
// }

// auto Interface::indexPhaseWithSpecies(Index ispecies) const -> Index
// {
//     Assert(ispecies < numSpecies(), "Cannot get the index of the phase with a given species.",
//         "The given species index `" + std::to_string(ispecies) + "` is out of range.");
//     const Index size = numPhases();
//     unsigned counter = 0;
//     for(unsigned i = 0; i < size; ++i)
//     {
//         counter += numSpeciesInPhase(i);
//         if(counter > ispecies) return i;
//     }
//     return size;
// }

// auto Interface::indexFirstSpeciesInPhase(Index iphase) const -> Index
// {
//     Assert(iphase < numPhases(), "Cannot get the index of first species in a given phase.",
//         "The given phase index `" + std::to_string(iphase) + "` is out of range.");
//     Index counter = 0;
//     for(unsigned i = 0; i < iphase; ++i)
//         counter += numSpeciesInPhase(i);
//     return counter;
// }

// auto Interface::system() const -> ChemicalSystem
// {
//     const unsigned nelements = numElements();
//     const unsigned nspecies = numSpecies();
//     const unsigned nphases = numPhases();

//     // Create a clone of the abstract Interface instance
//     // to be used in the lambda functions below.
//     std::shared_ptr<Interface> interface = clone();

//     // Create the Element instances
//     std::vector<Element> elements(nelements);
//     for(unsigned i = 0; i < nelements; ++i)
//     {
//         elements[i].setName(elementName(i));
//         elements[i].setMolarMass(elementMolarMass(i));
//     }

//     // Create the Species instances
//     std::vector<Species> species(nspecies);
//     for(unsigned i = 0; i < nspecies; ++i)
//     {
//         species[i].setName(speciesName(i));
//         species[i].setFormula(speciesName(i));
//         species[i].setElements(elementsInSpecies(*interface, elements, i));
//     }

//     // Create the Phase instances
//     std::vector<Phase> phases(nphases);
//     for(unsigned i = 0; i < nphases; ++i)
//     {
//         // Create the PhaseThermoModel function for the phase
//         PhaseThermoModel phase_thermo_model = [](PhaseThermoModelResult& res, real T, real P) -> void
//         {
//             RuntimeError("Could not evaluate the thermodynamic model of the phase.",
//                 "This phase was construted without a thermodynamic model.");
//         };

//         // Create the ActivityPropsFn function for the phase
//         ActivityPropsFn phase_chemical_model = [](ActivityPropsRef res, real T, real P, VectorXrConstRef n) -> void
//         {
//             RuntimeError("Could not evaluate the chemical model of phase.",
//                 "This phase was construted without a thermodynamic model.");
//         };

//         phases[i].setName(phaseName(i));
//         phases[i].setSpecies(speciesInPhase(*interface, species, i));
//         phases[i].setThermoModel(phase_thermo_model);
//         phases[i].setChemicalModel(phase_chemical_model);
//     }

//     // Create the ThermoModel function for the chemical system
//     ThermoModel thermo_model = [=](ThermoModelResult& res, real T, real P) -> void
//     {
//         interface->properties(res, T, P);
//     };

//     // Create the ChemicalModel function for the chemical system
//     ChemicalModel chemical_model = [=](ChemicalModelResult& res, real T, real P, VectorXrConstRef n) -> void
//     {
//         interface->properties(res, T, P, n);
//     };

//     // Create the ChemicalSystem instance
//     ChemicalSystem system(phases, thermo_model, chemical_model);

//     return system;
// }

// auto Interface::state(const ChemicalSystem& system) const -> ChemicalState
// {
//     ChemicalState state(system);
//     state.setTemperature(temperature());
//     state.setPressure(pressure());
//     VectorXr n = speciesAmounts();
//     n = (n.array() > 0.0).select(n, 1e-20);
//     state.setSpeciesAmounts(n);
//     return state;
// }

// Interface::operator ChemicalSystem() const
// {
//     return system();
// }

// } // namespace Reaktoro
