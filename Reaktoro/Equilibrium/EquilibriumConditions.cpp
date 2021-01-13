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

EquilibriumConditions::EquilibriumConditions(const EquilibriumSpecs& specs)
: mspecs(specs), mparamsmaster(specs.details().params)
{
}

auto EquilibriumConditions::temperature(real value, String unit) -> void
{
    error(!mparamsmaster.exists("T"), "Cannot set temperature for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "K");
    mparams.set("T", value);
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    error(!mparamsmaster.exists("P"), "Cannot set pressure for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "Pa");
    mparams.set("P", value);
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    error(!mparamsmaster.exists("V"), "Cannot set volume for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "m3");
    mparams.set("V", value);
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    error(!mparamsmaster.exists("U"), "Cannot set internal energy for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J");
    mparams.set("U", value);
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    error(!mparamsmaster.exists("H"), "Cannot set enthalpy for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J");
    mparams.set("H", value);
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    error(!mparamsmaster.exists("G"), "Cannot set Gibbs energy for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J");
    mparams.set("G", value);
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    error(!mparamsmaster.exists("A"), "Cannot set Helmholtz energy for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J");
    mparams.set("A", value);
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    error(!mparamsmaster.exists("S"), "Cannot set entropy for the equilibrium "
        "calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J/K");
    mparams.set("S", value);
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto paramname = "u[" + substance + "]";
    error(!mparamsmaster.exists(paramname), "Cannot set the chemical potential of ", substance, " for "
        "the equilibrium calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "J/mol");
    mparams.set(paramname, value);
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    error(!mparamsmaster.exists(paramname), "Cannot set the activity of ", species, " for "
        "the equilibrium calculation because it is not a registered input parameter in the specifications.");
    mparams.set(paramname, value);
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    error(!mparamsmaster.exists(paramname), "Cannot set the activity of ", species, " for "
        "the equilibrium calculation because it is not a registered input parameter in the specifications.");
    mparams.set(paramname, value * ln10);
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto paramname = "lnActivity[" + species + "]";
    error(!mparamsmaster.exists(paramname), "Cannot set the activity of ", species, " for "
        "the equilibrium calculation because it is not a registered input parameter in the specifications.");
    mparams.set(paramname, log(value));
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto paramname = "f[" + gas + "]";
    error(!mparamsmaster.exists(paramname), "Cannot set the fugacity of ", gas, " for "
        "the equilibrium calculation because it is not a registered input parameter in the specifications.");
    value = units::convert(value, unit, "bar");
    mparams.set(paramname, value);
}

auto EquilibriumConditions::pH(real value) -> void
{
    error(!mparamsmaster.exists("pH"), "Cannot set pH for the equilibrium calculation "
        "because it is not a registered input parameter in the specifications.");
    mparams.set("pH", value);
}

auto EquilibriumConditions::pMg(real value) -> void
{
    error(!mparamsmaster.exists("pMg"), "Cannot set pMg for the equilibrium calculation "
        "because it is not a registered input parameter in the specifications.");
    mparams.set("pMg", value);
}

auto EquilibriumConditions::pE(real value) -> void
{
    error(!mparamsmaster.exists("pE"), "Cannot set pE for the equilibrium calculation "
        "because it is not a registered input parameter in the specifications.");
    mparams.set("pE", value);
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    error(!mparamsmaster.exists("Eh"), "Cannot set Eh for the equilibrium calculation "
        "because it is not a registered input parameter in the specifications.");
    mparams.set("Eh", value);
}

auto EquilibriumConditions::system() const -> const ChemicalSystem&
{
    return mspecs.system();
}

auto EquilibriumConditions::specs() const -> const EquilibriumSpecs&
{
    return mspecs;
}

auto EquilibriumConditions::params() const -> const Params&
{
    return mparams;
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
