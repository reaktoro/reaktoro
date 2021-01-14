// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "EquilibriumConditions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>

namespace Reaktoro {
namespace {

/// Throw an error if an input parameter has not been registered in the equilibrium specifications.
/// @param params The list of registered input parameters in the equilibrium specifications (e.g., {"T", "P", "pH"}).
/// @param param The input parameter that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredParam(const Strings& params, const String& param, const String& propertymsg) -> void
{
    const auto registered = contains(params, param);
    error(!registered, "Cannot set ", propertymsg, " for the equilibrium calculation "
        "because it is not a registered input parameter in the equilibrium specifications.");
}

} // namespace

EquilibriumConditions::EquilibriumConditions(const EquilibriumSpecs& specs)
: m_system(specs.system()), m_parameters(specs.namesParameters())
{
}

auto EquilibriumConditions::temperature(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "T", "temperature");
    value = units::convert(value, unit, "K");
    m_params.set("T", value);
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "P", "pressure");
    value = units::convert(value, unit, "Pa");
    m_params.set("P", value);
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "V", "volume");
    value = units::convert(value, unit, "m3");
    m_params.set("V", value);
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "U", "internal energy");
    value = units::convert(value, unit, "J");
    m_params.set("U", value);
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "H", "enthalpy");
    value = units::convert(value, unit, "J");
    m_params.set("H", value);
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "G", "Gibbs energy");
    value = units::convert(value, unit, "J");
    m_params.set("G", value);
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "A", "Helmholtz energy");
    value = units::convert(value, unit, "J");
    m_params.set("A", value);
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "S", "entropy");
    value = units::convert(value, unit, "J/K");
    m_params.set("S", value);
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto paramname = "u[" + substance + "]";
    throwErrorIfNotRegisteredParam(m_parameters, paramname, "the chemical potential of " + substance);
    value = units::convert(value, unit, "J/mol");
    m_params.set(paramname, value);
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, paramname, "the activity of " + species);
    m_params.set(paramname, value);
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, paramname, "the activity of " + species);
    m_params.set(paramname, value * ln10);
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, paramname, "the activity of " + species);
    m_params.set(paramname, log(value));
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto paramname = "f[" + gas + "]";
    throwErrorIfNotRegisteredParam(m_parameters, paramname, "the fugacity of " + gas);
    value = units::convert(value, unit, "bar");
    m_params.set(paramname, value);
}

auto EquilibriumConditions::pH(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pH", "pH");
    m_params.set("pH", value);
}

auto EquilibriumConditions::pMg(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pMg", "pMg");
    m_params.set("pMg", value);
}

auto EquilibriumConditions::pE(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pE", "pE");
    m_params.set("pE", value);
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "Eh", "Eh");
    m_params.set("Eh", value);
}

auto EquilibriumConditions::initialComponentAmounts(VectorXrConstRef values) -> void
{
    const auto numcomponents = m_system.elements().size() + 1;
    error(values.size() != 0 && values.size() != numcomponents, "The number of conservative "
        "components is ", numcomponents, " but only ", values.size(), " values "
        "have been given for their amounts.");
    m_initial_component_amounts = values;
}

auto EquilibriumConditions::initialComponentAmounts() const -> VectorXrConstRef
{
    return m_initial_component_amounts;
}

auto EquilibriumConditions::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto EquilibriumConditions::params() const -> const Params&
{
    return m_params;
}


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

} // namespace Reaktoro
