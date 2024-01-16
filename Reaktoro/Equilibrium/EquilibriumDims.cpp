// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

EquilibriumDims::EquilibriumDims(EquilibriumSpecs const& specs)
{
    auto const& system = specs.system();

    Ne = system.elements().size();
    Nn = system.species().size();
    Np = specs.numControlVariablesP();
    Nq = specs.numControlVariablesQ();
    Nv = specs.numEquationConstraints();
    Nr = specs.numReactivityConstraints();
    Nb = 1 + Ne;
    Nc = 1 + Ne + Nr;
    Nt = specs.numTitrants();
    Nx = Nn + Nq;
    Nu = Nn + Np + Nq;
    Nw = specs.numInputs();

    errorif(Np != Nv,
        "The number of introduced p control variables (e.g., temperature, pressure, amounts of explicit titrants, custom variables) is ", Np, ". "
        "The number of introduced equation constraints is ", Nv, ". "
        "These two numbers must be equal, otherwise the chemical equilibrium problem cannot be solved. "
        "Modify your chemical equilibrium specifications so that this requirement is satisfied. ");
}

} // namespace Reaktoro
