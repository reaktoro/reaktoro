// Reaktor is a C++ library for computational reaction modelling.
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

#include "ReactionUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/ThermoProperty.hpp>
#include <Reaktor/Common/ThermoProperties.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/MultiphaseUtils.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

auto numSpecies(const Reaction& reaction) -> unsigned
{
	return reaction.species().size();
}

auto containsSpecies(const Reaction& reaction, const std::string& species) -> bool
{
	return indexSpecies(reaction, species) < numSpecies(reaction);
}

auto indexSpecies(const Reaction& reaction, const std::string& species) -> Index
{
	const auto& begin = reaction.species().begin();
	const auto& end = reaction.species().end();
	return std::find(begin, end, species) - begin;
}

auto indicesPhasesInReaction(const Multiphase& multiphase, const Reaction& reaction) -> Indices
{
	return indicesPhasesWithSpecies(multiphase, reaction.indices());
}

auto indicesReactionsWithSpecies(const Reactions& reactions, const Index& ispecies) -> Indices
{
	Indices indices;
	for(unsigned i = 0; i < reactions.size(); ++i)
		if(contained(ispecies, reactions[i].indices()))
			indices.push_back(i);
	return indices;
}

auto stoichiometry(const Reaction& reaction, const std::string& species) -> double
{
	const Index index = indexSpecies(reaction, species);
	return index < numSpecies(reaction) ? reaction.stoichiometries()[index] : 0.0;
}

//auto equilibriumConstant(const Multiphase& multiphase, const Reaction& reaction) -> EquilibriumConstant
//{
//	// The species in the chemical system
//    const std::vector<Species>& species = multiphase.species();
//
//    // The stoichiometries of the reacting species
//    const std::vector<double> stoichiometries = reaction.stoichiometries();
//
//    // The number of participating species in the reaction
//    const unsigned num_species = numSpecies(reaction);
//
//    // The universal gas constant (in units of J/(mol*K))
//    const double R = universalGasConstant;
//
//    // Collect the chemical potential functions of the reacting species
//    std::vector<ThermoPropertyFunction> mu;
//    for(Index i : reaction.indices())
//        mu.push_back(species[i].thermoModel().gibbs_energy);
//
//    // Define the equilibrium constant function
//    EquilibriumConstant kappa = [=](double T, double P)
//    {
//        double sum = 0.0;
//        for(unsigned i = 0; i < num_species; ++i)
//            sum += stoichiometries[i] * mu[i](T, P).val;
//        return std::exp(-sum/(R*T));
//    };
//
//    return kappa;
//}

template<typename PropertyFunction>
auto properties(const Reactions& reactions, double T, double P, PropertyFunction func) -> ThermoProperties
{
    const unsigned nreactions = reactions.size();
    Vector val(nreactions), ddt(nreactions), ddp(nreactions);
    for(unsigned i = 0; i < nreactions; ++i)
    {
        ThermoProperty prop = func(reactions[i], T, P);
        val[i] = prop.val;
        ddt[i] = prop.ddt;
        ddp[i] = prop.ddp;
    }
    return ThermoProperties(std::move(val), std::move(ddt), std::move(ddp));
}

auto equilibriumConstant(const Reaction& reaction, double T, double P) -> ThermoProperty
{
	return reaction.thermoModel().lnk(T, P);
}

auto equilibriumConstants(const Reactions& reactions, double T, double P) -> ThermoProperties
{
	return properties(reactions, T, P, equilibriumConstant);
}

auto rate(const Reaction& reaction, double T, double P, const Vector& n, const ThermoVector& a) -> ThermoScalar
{
	return reaction.kineticsModel().rate(T, P, n, a);
}

auto rates(const Reactions& reactions, double T, double P, const Vector& n, const ThermoVector& a) -> ThermoVector
{
	const unsigned nreactions = reactions.size();
	const unsigned nspecies = n.size();
	ThermoVector res(nreactions, nspecies);
	for(unsigned i = 0; i < nreactions; ++i)
		res.row(i) = rate(reactions[i], T, P, n, a);
	return res;
}

auto reactionQuotient(const Reaction& reaction, const ThermoVector& a) -> ThermoScalar
{
	const unsigned nspecies = a.val.size();

	ThermoScalar Q;
	Q.val = 1.0;
	Q.ddn = zeros(nspecies);

	const auto& stoichiometries = reaction.stoichiometries();
	const auto& indices = reaction.indices();

	for(unsigned i = 0; i < numSpecies(reaction); ++i)
	{
		const double vi = stoichiometries[i];
		const double ai = a.val[indices[i]];
		Q.val *= std::pow(ai, vi);
	}

	for(unsigned i = 0; i < numSpecies(reaction); ++i)
	{
		const double vi = stoichiometries[i];
		const double ai = a.val[indices[i]];
		Q.ddn += Q.val * vi/ai * a.ddn.row(indices[i]);
	}

	return Q;
}

auto reactionQuotients(const Reactions& reactions, const ThermoVector& a) -> ThermoVector
{
	const unsigned nreactions = reactions.size();
	const unsigned nspecies = a.val.size();
	ThermoVector res(nreactions, nspecies);
	for(unsigned i = 0; i < nreactions; ++i)
		res.row(i) = reactionQuotient(reactions[i], a);
	return res;
}

auto stoichiometricMatrix(const Multiphase& multiphase, const Reactions& reactions) -> Matrix
{
	const unsigned nrows = numSpecies(multiphase);
	const unsigned ncols = reactions.size();
	Matrix res(nrows, ncols);
	const auto& species = multiphase.species();
    for(unsigned i = 0; i < nrows; ++i)
        for(unsigned j = 0; j < ncols; ++j)
            res(i, j) = stoichiometry(reactions[i], species[j].name());
    return res;
}

} // namespace Reaktor
