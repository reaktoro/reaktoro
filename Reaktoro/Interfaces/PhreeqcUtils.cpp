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

#include "PhreeqcUtils.hpp"

namespace Reaktoro {
namespace PhreeqcUtils {

auto elements(const PhreeqcSpecies* species) -> std::map<PhreeqcElement*, double>
{
    std::map<element*, double> elements;
    for(auto iter = species->next_elt; iter->elt != nullptr; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto elements(const PhreeqcPhase* phase) -> std::map<PhreeqcElement*, double>
{
    std::map<element*, double> elements;
    for(auto iter = phase->next_elt; iter->elt != nullptr; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto stoichiometry(std::string element, const PhreeqcSpecies* species) -> double
{
    // Return species charge if element is electrical charge Z
    if(element == "Z") return species->z;

    // Iterate over all elements in the species
    for(auto iter = species->next_elt; iter->elt != nullptr; iter++)
        if(iter->elt->name == element)
            return iter->coef;
    return 0.0;
}

auto stoichiometry(std::string element, const PhreeqcPhase* phase) -> double
{
    // Iterate over all elements in the phase (gasesous or mineral species)
    for(auto iter = phase->next_elt; iter->elt != nullptr; iter++)
        if(iter->elt->name == element)
            return iter->coef;
    return 0.0;
}

auto isGaseousSpecies(const PhreeqcPhase* phase) -> bool
{
	std::string name(phase->name);
	return name.find("(g)") < name.size();
}

auto isMineralSpecies(const PhreeqcPhase* phase) -> bool
{
	return !isGaseousSpecies(phase);
}

auto reactionEquation(const PhreeqcSpecies* species) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;

    // Check if there is any reaction defined by this species
    if(species->rxn_x)
    {
        // Iterate over all reactants
        for(auto iter = species->rxn_x->token; iter->name != nullptr; iter++)
            equation.emplace(iter->name, iter->coef);

        // Check if the reaction is defined by only one species
        if(equation.size() <= 1)
            equation.clear();
    }

    return equation;
}

auto reactionEquation(const PhreeqcPhase* phase) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;
    for(auto iter = phase->rxn_x->token; iter->name != nullptr; iter++)
        equation.emplace(iter->name, iter->coef);
    return equation;
}

auto index(std::string name, const std::vector<PhreeqcSpecies*>& species) -> unsigned
{
    auto compare = [=](PhreeqcSpecies* entry) { return entry->name == name; };
    return std::find_if(species.begin(), species.end(), compare) - species.begin();
}

auto index(std::string name, const std::vector<PhreeqcPhase*>& phases) -> unsigned
{
    auto compare = [=](PhreeqcPhase* entry) { return entry->name == name; };
    return std::find_if(phases.begin(), phases.end(), compare) - phases.begin();
}

auto speciesAmounts(const std::vector<PhreeqcSpecies*>& species) -> Vector
{
    const unsigned size = species.size();
    Vector n(size);
    for(unsigned i = 0; i < size; ++i)
        n[i] = species[i]->moles;
    return n;
}

auto speciesAmounts(const std::vector<PhreeqcPhase*>& phases) -> Vector
{
    const unsigned size = phases.size();
    Vector n(size);
    for(unsigned i = 0; i < size; ++i)
        n[i] = phases[i]->moles_x;
    return n;
}

auto lnEquilibriumConstant(const double* l_logk, double T, double P) -> double
{
    // The universal gas constant (in units of kJ/(K*mol))
    const double R = 8.31470e-3;
    const double ln_10 = std::log(10.0);

    // Calculate log10(k) for this temperature and pressure
    double log10_k = l_logk[logK_T0]
        - l_logk[delta_h] * (298.15 - T)/(ln_10*R*T*298.15)
        + l_logk[T_A1]
        + l_logk[T_A2]*T
        + l_logk[T_A3]/T
        + l_logk[T_A4]*std::log10(T)
        + l_logk[T_A5]/(T*T)
        + l_logk[T_A6]*(T*T);

    return log10_k * ln_10;
}

auto lnEquilibriumConstant(const PhreeqcSpecies* species, double T, double P) -> double
{
    return lnEquilibriumConstant(species->rxn_x->logk, T, P);
}

auto lnEquilibriumConstant(const PhreeqcPhase* phase, double T, double P) -> double
{
    return lnEquilibriumConstant(phase->rxn_x->logk, T, P);
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
