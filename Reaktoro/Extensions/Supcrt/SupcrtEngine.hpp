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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoProps.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtParams.hpp>

namespace Reaktoro {

/// The class used for standard thermodynamic property calculations based on SUPCRT.
/// @ingroup SupcrtExtension
class SupcrtEngine
{
public:
    /// Construct a default SupcrtEngine object.
    SupcrtEngine();

    /// Return the standard thermodynamic properties of solvent water using the HKF model.
    auto props(real T, real P, const SupcrtParamsAqueousSolventHKF& params) const -> StandardThermoProps;

    /// Return the standard thermodynamic properties of an aqueous solute using the HKF model.
    auto props(real T, real P, const SupcrtParamsAqueousSoluteHKF& params) const -> StandardThermoProps;

    /// Return the standard thermodynamic properties of a fluid species using the Maier-Kelly model.
    auto props(real T, real P, const SupcrtParamsMaierKelly& params) const -> StandardThermoProps;

    /// Return the standard thermodynamic properties of a mineral species using the Maier-Kelly-HKF model.
    auto props(real T, real P, const SupcrtParamsMaierKellyHKF& params) const -> StandardThermoProps;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
