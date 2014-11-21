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
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/MultiphaseUtils.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {
//
//auto numSpecies(const Reaction& reaction) -> unsigned
//{
//	return reaction.species().size();
//}
//
//auto containsSpecies(const Reaction& reaction, const std::string& species) -> bool
//{
//    return speciesIndex(reaction, species) < numSpecies(reaction);
//}
//
//auto speciesIndex(const Reaction& reaction, const std::string& species) -> Index
//{
//    return find(species, reaction.species());
//}
//
//auto phaseIndicesInReaction(const Multiphase& multiphase, const Reaction& reaction) -> Indices
//{
//	return phaseIndicesWithSpecies(multiphase, reaction.indices());
//}
//
//auto indicesReactionsWithSpecies(const Reactions& reactions, const Index& ispecies) -> Indices
//{
//	Indices indices;
//	for(unsigned i = 0; i < reactions.size(); ++i)
//		if(contained(ispecies, reactions[i].indices()))
//			indices.push_back(i);
//	return indices;
//}
//
//auto stoichiometry(const Reaction& reaction, const std::string& species) -> double
//{
//	const Index index = speciesIndex(reaction, species);
//	return index < numSpecies(reaction) ? reaction.stoichiometries()[index] : 0.0;
//}
//
//auto auxThermoScalarFunction(const Reaction& reaction, const std::vector<ThermoScalarFunction>& f) -> ThermoScalarFunction
//{
//    // The stoichiometries of the reacting species
//    const std::vector<double> stoichiometries = reaction.stoichiometries();
//
//    // Define the thermo property function
//    ThermoScalarFunction thermofn = [=](double T, double P)
//    {
//        double res_val = 0.0;
//        double res_ddt = 0.0;
//        double res_ddp = 0.0;
//        for(unsigned i = 0; i < stoichiometries.size(); ++i)
//        {
//            ThermoScalar res = f[i](T, P);
//            res_val += stoichiometries[i] * res.val();
//            res_ddt += stoichiometries[i] * res.ddt();
//            res_ddp += stoichiometries[i] * res.ddp();
//        }
//        return ThermoScalar(res_val, res_ddt, res_ddp);
//    };
//
//    return thermofn;
//}
//
//auto thermoModel(const Multiphase& multiphase, const Reaction& reaction) -> ReactionThermoModel
//{
//    // The species in the chemical system
//    const std::vector<Species>& species = multiphase.species();
//
//    // The universal gas constant (in units of J/(mol*K))
//    const double R = universalGasConstant;
//
//    // Collect the thermodynamic property functions of all reactants
//    std::vector<ThermoScalarFunction> gibbs_fns;
//    std::vector<ThermoScalarFunction> enthalpy_fns;
//    std::vector<ThermoScalarFunction> entropy_fns;
//    std::vector<ThermoScalarFunction> ienergy_fns;
//    std::vector<ThermoScalarFunction> helmholtz_fns;
//
//    for(Index i : reaction.indices())
//    {
//        gibbs_fns.push_back(species[i].thermoModel().gibbs_energy);
//        enthalpy_fns.push_back(species[i].thermoModel().enthalpy);
//        entropy_fns.push_back(species[i].thermoModel().entropy);
//        ienergy_fns.push_back(species[i].thermoModel().internal_energy);
//        helmholtz_fns.push_back(species[i].thermoModel().helmholtz_energy);
//    }
//
//    // Set the thermodynamic property functions of the reaction
//    ReactionThermoModel thermo_model;
//    thermo_model.gibbs_energy     = auxThermoScalarFunction(reaction, gibbs_fns);
//    thermo_model.enthalpy         = auxThermoScalarFunction(reaction, enthalpy_fns);
//    thermo_model.entropy          = auxThermoScalarFunction(reaction, entropy_fns);
//    thermo_model.internal_energy  = auxThermoScalarFunction(reaction, ienergy_fns);
//    thermo_model.helmholtz_energy = auxThermoScalarFunction(reaction, helmholtz_fns);
//
//    // Set the equilibrium constant thermodynamic property of the reaction
//    ThermoScalarFunction gibbsfn = thermo_model.gibbs_energy;
//
//    thermo_model.lnk = [=](double T, double P)
//    {
//        ThermoScalar g = gibbsfn(T, P);
//        const double lnk_val = -g.val()/(R*T);
//        const double lnk_ddt = -g.ddt()/(R*T) + g.val()/(R*T*T);
//        const double lnk_ddp = -g.ddp()/(R*T);
//
//        return ThermoScalar(lnk_val, lnk_ddt, lnk_ddp);
//    };
//
//    return thermo_model;
//}
//
//template<typename PropertyFunction>
//auto properties(const Reactions& reactions, double T, double P, PropertyFunction func) -> ThermoVector
//{
//    const unsigned nreactions = reactions.size();
//    ThermoVector res(nreactions);
//    for(unsigned i = 0; i < nreactions; ++i)
//        res.row(i) = func(reactions[i], T, P);
//    return res;
//}
//
//auto equilibriumConstant(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//	return reaction.thermoModel().lnk(T, P);
//}
//
//auto equilibriumConstants(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//	return properties(reactions, T, P, equilibriumConstant);
//}
//
//auto entropy(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//    return reaction.thermoModel().entropy(T, P);
//}
//
//auto entropies(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//    return properties(reactions, T, P, entropy);
//}
//
//auto helmholtzEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//    return reaction.thermoModel().helmholtz_energy(T, P);
//}
//
//auto helmholtzEnergies(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//    return properties(reactions, T, P, helmholtzEnergy);
//}
//
//auto internalEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//    return reaction.thermoModel().internal_energy(T, P);
//}
//
//auto internalEnergies(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//    return properties(reactions, T, P, internalEnergy);
//}
//
//auto enthalpy(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//    return reaction.thermoModel().enthalpy(T, P);
//}
//
//auto enthalpies(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//    return properties(reactions, T, P, enthalpy);
//}
//
//auto gibbsEnergy(const Reaction& reaction, double T, double P) -> ThermoScalar
//{
//    return reaction.thermoModel().gibbs_energy(T, P);
//}
//
//auto gibbsEnergies(const Reactions& reactions, double T, double P) -> ThermoVector
//{
//    return properties(reactions, T, P, gibbsEnergy);
//}
//
//auto rate(const Reaction& reaction, double T, double P, const Vector& n, const ChemicalVector& a) -> ChemicalScalar
//{
//	return reaction.kineticsModel().rate(T, P, n, a);
//}
//
//auto rates(const Reactions& reactions, double T, double P, const Vector& n, const ChemicalVector& a) -> ChemicalVector
//{
//	const unsigned nreactions = reactions.size();
//	const unsigned nspecies = n.size();
//	ChemicalVector r(nreactions, nspecies);
//	for(unsigned i = 0; i < nreactions; ++i)
//	    r.row(i) = rate(reactions[i], T, P, n, a);
//	return r;
//}
//
//auto reactionQuotient(const Reaction& reaction, const ChemicalVector& a) -> ChemicalScalar
//{
//	const unsigned nspecies = a.val().size();
//
//	double Q_val = 1.0;
//	Vector Q_ddn = zeros(nspecies);
//
//	const auto& stoichiometries = reaction.stoichiometries();
//	const auto& indices = reaction.indices();
//
//	for(unsigned i = 0; i < numSpecies(reaction); ++i)
//	{
//		const double vi = stoichiometries[i];
//		const double ai = a.val()[indices[i]];
//		Q_val *= std::pow(ai, vi);
//	}
//
//	for(unsigned i = 0; i < numSpecies(reaction); ++i)
//	{
//		const double vi = stoichiometries[i];
//		const double ai = a.val()[indices[i]];
//		Q_ddn += Q_val * vi/ai * a.ddn().row(indices[i]);
//	}
//
//	return {Q_val, 0.0, 0.0, Q_ddn};
//}
//
//auto reactionQuotients(const Reactions& reactions, const ChemicalVector& a) -> ChemicalVector
//{
//	const unsigned nreactions = reactions.size();
//	const unsigned nspecies = a.val().size();
//	ChemicalVector res(nreactions, nspecies);
//	for(unsigned i = 0; i < nreactions; ++i)
//		res.row(i) = reactionQuotient(reactions[i], a);
//	return res;
//}
//
//auto stoichiometricMatrix(const Multiphase& multiphase, const Reactions& reactions) -> Matrix
//{
//	const unsigned nrows = numSpecies(multiphase);
//	const unsigned ncols = reactions.size();
//	Matrix res(nrows, ncols);
//	const auto& species = multiphase.species();
//    for(unsigned i = 0; i < nrows; ++i)
//        for(unsigned j = 0; j < ncols; ++j)
//            res(i, j) = stoichiometry(reactions[i], species[j].name());
//    return res;
//}

} // namespace Reaktor
