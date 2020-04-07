// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "ActivityModelIdealAqueousSolution.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

auto ActivityModelIdealAqueousSolution::create(const SpeciesList& species) -> ActivityModelFn
{
    const auto iH2O = species.indexWithFormula("H2O");
    const auto MH2O = waterMolarMass;

    ActivityModelFn fn = [=](ActivityProps& props, real T, real P, ArrayXrConstRef x)
    {
        using std::log;
        const auto m = x/(MH2O * x[iH2O]); // molalities
        props.Vex  = 0.0;
        props.VexT = 0.0;
        props.VexP = 0.0;
        props.Gex  = 0.0;
        props.Hex  = 0.0;
        props.Cpex = 0.0;
        props.ln_g = 0.0;
        props.ln_a = m.log();
        props.ln_a[iH2O] = log(x[iH2O]);
    };

    return fn;
}

} // namespace Reaktoro
