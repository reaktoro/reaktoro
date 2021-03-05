// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

/// Throw an error if an input parameter has not been registered in the equilibrium specifications.
/// @param params The list of registered input parameters in the equilibrium specifications (e.g., {"T", "P", "pH"}).
/// @param pid The id of the parameter that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredParam(const Params& params, const String& pid, const String& propertymsg) -> void
{
    const auto registered = containsfn(params, RKT_LAMBDA(x, x.id() == pid));
    error(!registered, "Cannot set ", propertymsg, " for the equilibrium calculation "
        "because it is not a registered input parameter in the equilibrium specifications.");
}

} // namespace

EquilibriumConditions::EquilibriumConditions(const EquilibriumSpecs& specs)
: m_system(specs.system()), m_params(specs.params())
{}

auto EquilibriumConditions::temperature(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "T", "temperature");
    value = units::convert(value, unit, "K");
    m_params.get("T") = value;
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "P", "pressure");
    value = units::convert(value, unit, "Pa");
    m_params.get("P") = value;
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "V", "volume");
    value = units::convert(value, unit, "m3");
    m_params.get("V") = value;
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "U", "internal energy");
    value = units::convert(value, unit, "J");
    m_params.get("U") = value;
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "H", "enthalpy");
    value = units::convert(value, unit, "J");
    m_params.get("H") = value;
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "G", "Gibbs energy");
    value = units::convert(value, unit, "J");
    m_params.get("G") = value;
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "A", "Helmholtz energy");
    value = units::convert(value, unit, "J");
    m_params.get("A") = value;
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "S", "entropy");
    value = units::convert(value, unit, "J/K");
    m_params.get("S") = value;
}

auto EquilibriumConditions::startWith(String species, real value, String unit) -> void
{
    const auto ispecies = m_system.species().index(species);
    startWith(ispecies, value, unit);
}

auto EquilibriumConditions::startWith(Index ispecies, real value, String unit) -> void
{
    const auto size = m_system.species().size();
    errorif(ispecies >= size, "Given species index (", ispecies, ") is out-of-bounds (number of species is ", size, ").");
    const auto amount = detail::computeSpeciesAmount(m_system, ispecies, value, unit);
    m_initial_species_amounts.resize(size);
    m_initial_species_amounts[ispecies] = amount;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumConditions::startWith(const ChemicalState& state) -> void
{
    const auto n = state.speciesAmounts();
    const auto size = m_system.species().size();
    errorif(n.size() != size, "Given chemical state does not have compatible number of species (", n.size(), ") with number of species in the system (", size, ")");
    m_initial_species_amounts = n;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumConditions::startWithComponentAmounts(ArrayXrConstRef b) -> void
{
    const auto Nb = m_system.elements().size() + 1;
    errorif(b.size() != Nb, "The number of conservative components is ", Nb, " but only ", b.size(), " values have been given for their initial amounts.");
    m_initial_species_amounts.resize(0);
    m_initial_component_amounts = b;
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto pid = "u[" + substance + "]";
    throwErrorIfNotRegisteredParam(m_params, pid, "the chemical potential of " + substance);
    value = units::convert(value, unit, "J/mol");
    m_params.get(pid) = value;
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_params, pid, "the activity of " + species);
    m_params.get(pid) = value;
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_params, pid, "the activity of " + species);
    m_params.get(pid) = value * ln10;
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredParam(m_params, pid, "the activity of " + species);
    m_params.get(pid) = log(value);
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto pid = "f[" + gas + "]";
    throwErrorIfNotRegisteredParam(m_params, pid, "the fugacity of " + gas);
    value = units::convert(value, unit, "bar");
    m_params.get(pid) = value;
}

auto EquilibriumConditions::pH(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "pH", "pH");
    m_params.get("pH") = value;
}

auto EquilibriumConditions::pMg(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "pMg", "pMg");
    m_params.get("pMg") = value;
}

auto EquilibriumConditions::pE(real value) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "pE", "pE");
    m_params.get("pE") = value;
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    throwErrorIfNotRegisteredParam(m_params, "Eh", "Eh");
    m_params.get("Eh") = value;
}

auto EquilibriumConditions::initialSpeciesAmounts() const -> ArrayXrConstRef
{
    return m_initial_species_amounts;
}

auto EquilibriumConditions::initialComponentAmounts() const -> ArrayXr
{
    if(m_initial_component_amounts.size())
        return m_initial_component_amounts;

    errorif(m_initial_species_amounts.size() == 0, "Initial conditions for species amounts or components amounts have not been given.");

    const auto Wn = m_system.formulaMatrix();
    const auto n0 = m_initial_species_amounts.matrix();
    const auto b = Wn * n0;
    return b;
}

auto EquilibriumConditions::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto EquilibriumConditions::params() const -> Params
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
