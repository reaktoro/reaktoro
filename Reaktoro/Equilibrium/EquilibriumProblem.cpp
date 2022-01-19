// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "EquilibriumProblem.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace detail {

/// Return an EquilibriumSpecs object that represent the specifications of a Gibbs energy minimization problem.
inline auto createEquilibriumSpecs(const ChemicalSystem& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    return specs;
}

} // namespace detail

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: EquilibriumProblem(detail::createEquilibriumSpecs(system))
{}

EquilibriumProblem::EquilibriumProblem(const EquilibriumSpecs& specs)
: EquilibriumConditions(specs), EquilibriumRestrictions(specs.system())
{}

auto EquilibriumProblem::startWithTemperature(real value) -> void
{
    errorif(value <= 0.0, "EquilibriumProblem::startWithTemperature requires a positive temperature value in K, but the given value was ", value, " K.");
    m_initial_temperature = value;
}

auto EquilibriumProblem::startWithTemperature(real value, String unit) -> void
{
    auto converted = units::convert(value, unit, "K");
    errorif(converted <= 0.0, "EquilibriumProblem::startWithTemperature requires a positive temperature value in K, but the given value was ", value, " ", unit);
    m_initial_temperature = converted;
}

auto EquilibriumProblem::startWithPressure(real value) -> void
{
    errorif(value <= 0.0, "EquilibriumProblem::startWithPressure requires a positive pressure value in Pa, but the given value was ", value, " Pa.");
    m_initial_pressure = value;
}

auto EquilibriumProblem::startWithPressure(real value, String unit) -> void
{
    auto converted = units::convert(value, unit, "K");
    errorif(converted <= 0.0, "EquilibriumProblem::startWithPressure requires a positive pressure value in Pa, but the given value was ", value, " ", unit);
    m_initial_pressure = converted;
}

auto EquilibriumProblem::startWith(String species, real value, String unit) -> void
{
    const auto ispecies = system().species().index(species);
    startWith(ispecies, value, unit);
}

auto EquilibriumProblem::startWith(Index ispecies, real value, String unit) -> void
{
    const auto size = system().species().size();
    errorif(ispecies >= size, "EquilibriumProblem::startWith requires a valid species index (got index ", ispecies, ", which is greater or equal than the number of species", size, ")");
    const auto amount = detail::computeSpeciesAmount(system(), ispecies, value, unit);
    m_initial_species_amounts.resize(size);
    m_initial_species_amounts[ispecies] = amount;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumProblem::startWithSpeciesAmounts(ArrayXrConstRef n) -> void
{
    const auto size = system().species().size();
    errorif(n.size() != size, "EquilibriumProblem::startWithSpeciesAmounts requires an array with as many entries as there are chemical species in the system (expected size is ", size, ", but given was", n.size(), ")");
    m_initial_species_amounts = n;
    m_initial_component_amounts.resize(0);
}

auto EquilibriumProblem::startWithState(const ChemicalState& state) -> void
{
    startWithTemperature(state.temperature());
    startWithPressure(state.pressure());
    startWithSpeciesAmounts(state.speciesAmounts());
}

auto EquilibriumProblem::startWithComponentAmounts(ArrayXrConstRef b) -> void
{
    const auto Nb = system().elements().size() + 1;
    errorif(b.size() != Nb, "EquilibriumProblem::startWithComponentAmounts requires an array with as many entries as there are convervative components in the equilibrium problem (expected size is ", Nb, ", but given was", b.size(), ")");
    m_initial_component_amounts = b;
}

auto EquilibriumProblem::initialTemperature() const -> real
{
    return m_initial_temperature;
}

auto EquilibriumProblem::initialPressure() const -> real
{
    return m_initial_pressure;
}

auto EquilibriumProblem::initialSpeciesAmounts() const -> ArrayXrConstRef
{
    return m_initial_species_amounts;
}

auto EquilibriumProblem::initialComponentAmounts() const -> ArrayXr
{
    if(m_initial_component_amounts.size())
        return m_initial_component_amounts;

    errorif(m_initial_species_amounts.size() == 0,
        "While executing EquilibriumProblem::initialComponentAmounts, it was found "
        "that initial conditions for species or component amounts have not been given. "
        "You need to use one of the methods below (check their overloaded versions):\n"
        " * EquilibriumProblem::startWith\n"
        " * EquilibriumProblem::startWithComponentAmounts");

    const auto Wn = system().formulaMatrix();
    const auto n0 = m_initial_species_amounts.matrix();
    const auto b = Wn * n0;
    return b;
}
} // namespace Reaktoro
