//// Reaktoro is a unified framework for modeling chemically reactive systems.
////
//// Copyright (C) 2014-2015 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#pragma once
//
//// Reaktoro includes
//#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>
//
//namespace Reaktoro {
//
//// Forward declarations
//class MineralMixture;
//
///// Return an equation of state for a binary mineral solid solution based on Van Laar model.
///// The Van Laar model calculates the activity coefficient of the end-members in a
///// solid solution using the equations:
///// @todo Write here the equations todo
///// @param mixture The mineral mixture
///// @param a The size parameters for the solid solution end-members
///// @param W The binary interaction parameters for the solid solution end-members
///// @return The equation of state function for the mineral phase
///// @see MineralMixture, MineralChemicalModel
//auto mineralChemicalModelVanLaar(const MineralMixture& mixture, const Vector& a, const Matrix& W) -> PhaseChemicalModel;
//
//} // namespace Reaktoro
