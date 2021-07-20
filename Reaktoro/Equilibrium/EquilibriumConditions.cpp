// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

/// Throw an error if an input variable has not been registered in the equilibrium specifications.
/// @param inputs The list of registered input variables in the equilibrium specifications (e.g., {"T", "P", "pH"}).
/// @param wid The id of the input variable that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredInput(const Strings& inputs, const String& wid, const String& propertymsg) -> void
{
    const auto registered = contains(inputs, wid);
    errorif(!registered, "Cannot set ", propertymsg, " for the equilibrium calculation "
        "because it is not a registered input variable in the equilibrium specifications.");
}

} // namespace

EquilibriumConditions::EquilibriumConditions(const EquilibriumSpecs& specs)
: m_system(specs.system()), m_inputs(specs.inputs())
{
    // Initialize the values of the input variables to zero
    m_inputs_values = zeros(specs.numInputs());

    // Initialize the values of the input variables that are model parameters to their current values
    m_inputs_values(specs.indicesParams()) = VectorXr(specs.params());
}

auto EquilibriumConditions::temperature(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "T", "temperature");
    value = units::convert(value, unit, "K");
    const auto idx = index(m_inputs, "T");
    m_inputs_values[idx] = value;
    startWithTemperature(value);
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "P", "pressure");
    value = units::convert(value, unit, "Pa");
    const auto idx = index(m_inputs, "P");
    m_inputs_values[idx] = value;
    startWithPressure(value);
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "V", "volume");
    value = units::convert(value, unit, "m3");
    const auto idx = index(m_inputs, "V");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "U", "internal energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "U");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "H", "enthalpy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "H");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "G", "Gibbs energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "G");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "A", "Helmholtz energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "A");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "S", "entropy");
    value = units::convert(value, unit, "J/K");
    const auto idx = index(m_inputs, "S");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::startWithTemperature(real value) -> void
{
    errorif(value <= 0.0, "EquilibriumConditions::startWithTemperature requires a positive temperature value in K, but the given value was ", value, " K.");
    m_initial_temperature = value;
}

auto EquilibriumConditions::startWithTemperature(real value, String unit) -> void
{
    auto converted = units::convert(value, unit, "K");
    errorif(converted <= 0.0, "EquilibriumConditions::startWithTemperature requires a positive temperature value in K, but the given value was ", value, " ", unit);
    m_initial_temperature = converted;
}

auto EquilibriumConditions::startWithPressure(real value) -> void
{
    errorif(value <= 0.0, "EquilibriumConditions::startWithPressure requires a positive pressure value in Pa, but the given value was ", value, " Pa.");
    m_initial_pressure = value;
}

auto EquilibriumConditions::startWithPressure(real value, String unit) -> void
{
    auto converted = units::convert(value, unit, "K");
    errorif(converted <= 0.0, "EquilibriumConditions::startWithPressure requires a positive pressure value in Pa, but the given value was ", value, " ", unit);
    m_initial_pressure = converted;
}

auto EquilibriumConditions::startWith(String species, real value, String unit) -> void
{
    const auto ispecies = m_system.species().index(species);
    startWith(ispecies, value, unit);
}

auto EquilibriumConditions::startWith(Index ispecies, real value, String unit) -> void
{
    const auto size = m_system.species().size();
    errorif(ispecies >= size, "EquilibriumConditions::startWith requires a valid species index (got index ", ispecies, ", which is greater or equal than the number of species", size, ")");
    const auto amount = detail::computeSpeciesAmount(m_system, ispecies, value, unit);
    m_initial_species_amounts.resize(size);
    m_initial_species_amounts[ispecies] = amount;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumConditions::startWithSpeciesAmounts(ArrayXrConstRef n) -> void
{
    const auto size = m_system.species().size();
    errorif(n.size() != size, "EquilibriumConditions::startWithSpeciesAmounts requires an array with as many entries as there are chemical species in the system (expected size is ", size, ", but given was", n.size(), ")");
    m_initial_species_amounts = n;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumConditions::startWithState(const ChemicalState& state) -> void
{
    startWithTemperature(state.temperature());
    startWithPressure(state.pressure());
    startWithSpeciesAmounts(state.speciesAmounts());
}

auto EquilibriumConditions::startWithComponentAmounts(ArrayXrConstRef b) -> void
{
    const auto Nb = m_system.elements().size() + 1;
    errorif(b.size() != Nb, "EquilibriumConditions::startWithComponentAmounts requires an array with as many entries as there are convervative components in the equilibrium problem (expected size is ", Nb, ", but given was", b.size(), ")");
    m_initial_component_amounts = b;
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto pid = "u[" + substance + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the chemical potential of " + substance);
    value = units::convert(value, unit, "J/mol");
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value * ln10;
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto pid = "lnActivity[" + species + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = log(value);
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto pid = "f[" + gas + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the fugacity of " + gas);
    value = units::convert(value, unit, "bar");
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pH(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pH", "pH");
    const auto idx = index(m_inputs, "pH");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pMg(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pMg", "pMg");
    const auto idx = index(m_inputs, "pMg");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pE(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pE", "pE");
    const auto idx = index(m_inputs, "pE");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "Eh", "Eh");
    const auto idx = index(m_inputs, "Eh");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::set(const String& input, const real& val) -> void
{
    const auto idx = index(m_inputs, input);
    const auto size = m_inputs.size();
    errorif(idx >= size, "There is no input variable with name `", input, "` in this EquilibriumConditions object.");
    m_inputs_values[idx] = val;
}

auto EquilibriumConditions::initialTemperature() const -> real
{
    return m_initial_temperature;
}

auto EquilibriumConditions::initialPressure() const -> real
{
    return m_initial_pressure;
}

auto EquilibriumConditions::initialSpeciesAmounts() const -> ArrayXrConstRef
{
    return m_initial_species_amounts;
}

auto EquilibriumConditions::initialComponentAmounts() const -> ArrayXr
{
    if(m_initial_component_amounts.size())
        return m_initial_component_amounts;

    errorif(m_initial_species_amounts.size() == 0,
        "While executing EquilibriumConditions::initialComponentAmounts, it was found "
        "that initial conditions for species or component amounts have not been given. "
        "You need to use one of the methods below (check their overloaded versions):\n"
        " * EquilibriumConditions::startWith\n"
        " * EquilibriumConditions::startWithComponentAmounts");

    const auto Wn = m_system.formulaMatrix();
    const auto n0 = m_initial_species_amounts.matrix();
    const auto b = Wn * n0;
    return b;
}

auto EquilibriumConditions::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto EquilibriumConditions::inputNames() const -> const Strings&
{
    return m_inputs;
}

auto EquilibriumConditions::inputValues() const -> VectorXrConstRef
{
    return m_inputs_values;
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
