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
//
//#pragma once
//
//// Reaktoro includes
//#include <Reaktoro/Core/ActivityModel.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class GeneralMixture;
//
///// Return an equation of state for a binary mineral solid solution based on Van Laar model.
///// The Van Laar model calculates the activity coefficient of the end-members in a
///// solid solution using the equations:
///// @todo Write here the equations todo
///// @param mixture The mineral mixture
///// @param a The size parameters for the solid solution end-members
///// @param W The binary interaction parameters for the solid solution end-members
///// @return The equation of state function for the mineral phase
///// @see GeneralMixture, MineralChemicalModel
//auto mineralChemicalModelVanLaar(const GeneralMixture& mixture, VectorXrConstRef a, MatrixXdConstRef W)-> ActivityModelFn;
//
//} // namespace Reaktoro
