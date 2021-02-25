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

#include "EquilibriumDims.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

EquilibriumDims::EquilibriumDims(const EquilibriumSpecs& specs)
{
    const auto& system = specs.system();

    Ne = system.elements().size() + 1;
    Nb = Ne; // TODO: Currently, this is chemical elements + electric charge. But we should change this when using EquilibriumReactions, where we will define the components (possibly fictitious ones if reactions are prevented in the equilibrium calculation).
    Nn = system.species().size();
    Np = specs.numControlVariables() - specs.numTitrantsImplicit();
    Nq = specs.numTitrantsImplicit();
    Nt = specs.numTitrants();
    Nx = Nn + Nq;
    Nu = Nn + Np + Nq;

    error(Np + Nq != specs.numConstraints(),
        "The number of introduced control variables (e.g., temperature, pressure, amounts of titrants) is ", Np + Nq, ". "
        "The number of introduced constraints (e.g., equation constraints and chemical potential constraints) is ", specs.numConstraints(), ". "
        "These two numbers must be equal, otherwise the chemical equilibrium problem cannot be solved. "
        "Modify your chemical equilibrium specifications so that this requirement is satisfied. ");
}

} // namespace Reaktoro
