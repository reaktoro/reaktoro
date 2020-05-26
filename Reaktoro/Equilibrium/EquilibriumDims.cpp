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

#include "EquilibriumDims.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>

namespace Reaktoro {

EquilibriumDims::EquilibriumDims(const EquilibriumConstraints& constraints)
{
    const auto& system = constraints.system();

    Ne  = system.elements().size() + 1;
    Nn  = system.species().size();
    Npe = constraints.data().econstraints.size();
    Npp = constraints.data().pconstraints.size();
    Np  = Npe + Npp;
    Nq  = constraints.data().uconstraints.size();
    Nir = constraints.data().restrictions.reactions_cannot_react.size();
    Ncv = constraints.data().controls.size() + Nq;
    Nx  = Nn + Np + Nq;
    Nc  = Ne + Nir;

    error(Np + Nq != Ncv,
        "The number of introduced control variables (", Ncv - Nq, ") "
        "using method EquilibriumConstraints::control must be equal to the "
        "number of functional constraints introduced with methods "
        "EquilibriumConstraints::until (", Npe, ") and "
        "EquilibriumConstraints::preserve (", Npp, ").");
}

} // namespace Reaktoro
