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
/// @param pid The id of the parameter that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredParam(const Strings& params, const String& pid, const String& propertymsg) -> void
{
    const auto registered = contains(params, pid);
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
    m_params.append("T", value);
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "P", "pressure");
    value = units::convert(value, unit, "Pa");
    m_params.append("P", value);
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "V", "volume");
    value = units::convert(value, unit, "m3");
    m_params.append("V", value);
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "U", "internal energy");
    value = units::convert(value, unit, "J");
    m_params.append("U", value);
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "H", "enthalpy");
    value = units::convert(value, unit, "J");
    m_params.append("H", value);
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "G", "Gibbs energy");
    value = units::convert(value, unit, "J");
    m_params.append("G", value);
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "A", "Helmholtz energy");
    value = units::convert(value, unit, "J");
    m_params.append("A", value);
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "S", "entropy");
    value = units::convert(value, unit, "J/K");
    m_params.append("S", value);
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto pid = "u[" + substance + "]";
    throwErrorIfNotRegisteredParam(m_parameters, pid, "the chemical potential of " + substance);
    value = units::convert(value, unit, "J/mol");
    m_params.append(pid, value);
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, pid, "the activity of " + species);
    m_params.append(pid, value);
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, pid, "the activity of " + species);
    m_params.append(pid, value * ln10);
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_parameters, pid, "the activity of " + species);
    m_params.append(pid, log(value));
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto pid = "f[" + gas + "]";
    throwErrorIfNotRegisteredParam(m_parameters, pid, "the fugacity of " + gas);
    value = units::convert(value, unit, "bar");
    m_params.append(pid, value);
}

auto EquilibriumConditions::pH(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pH", "pH");
    m_params.append("pH", value);
}

auto EquilibriumConditions::pMg(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pMg", "pMg");
    m_params.append("pMg", value);
}

auto EquilibriumConditions::pE(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "pE", "pE");
    m_params.append("pE", value);
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_parameters, "Eh", "Eh");
    m_params.append("Eh", value);
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
