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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Models/ActivityModels/Support/CubicEOS.hpp>

namespace Reaktoro {

/// The type for functions that construct a CubicEOS::BipModel for a fluid phase once its species are known.
/// @param specieslist The species in the fluid phase.
using CubicBipModelGenerator = Fn<CubicEOS::BipModel(SpeciesList const& specieslist)>;

/// Return the activity model for fluid phases based on the Van der Waals cubic equation of state.
auto ActivityModelVanDerWaals(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Redlich-Kwong cubic equation of state.
auto ActivityModelRedlichKwong(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Soave-Redlich-Kwong cubic equation of state.
auto ActivityModelSoaveRedlichKwong(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson (1978) cubic equation of state.
auto ActivityModelPengRobinson(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson (1976) cubic equation of state.
auto ActivityModelPengRobinson76(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson (1978) cubic equation of state.
auto ActivityModelPengRobinson78(CubicBipModelGenerator cbipmodel = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson (1976) with the binary interaction parameter model used in PHREEQC.
auto ActivityModelPengRobinsonPHREEQC() -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson (1978) with the binary interaction parameter model of Søreide and Whitson (1992).
auto ActivityModelPengRobinsonSoreideWhitson() -> ActivityModelGenerator;

/// Return the binary interaction parameter model for Peng-Robinson EOS (1976) equivalent to that used in PHREEQC.
auto CubicBipModelPHREEQC() -> CubicBipModelGenerator;

/// Return the binary interaction parameter model for Peng-Robinson EOS (1978) equivalent to that reported in Søreide and Whitson (1992).
auto CubicBipModelSoreideWhitson() -> CubicBipModelGenerator;

} // namespace Reaktoro
