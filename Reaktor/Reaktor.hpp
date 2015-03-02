// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <Reaktor/Activity/ActivityUtils.hpp>
#include <Reaktor/Activity/AqueousActivity.hpp>
#include <Reaktor/Activity/AqueousActivityDrummond.hpp>
#include <Reaktor/Activity/AqueousActivityDuanSun.hpp>
#include <Reaktor/Activity/AqueousActivityHKF.hpp>
#include <Reaktor/Activity/AqueousActivityIdeal.hpp>
#include <Reaktor/Activity/AqueousActivityPitzer.hpp>
#include <Reaktor/Activity/AqueousActivityRumpf.hpp>
#include <Reaktor/Activity/AqueousActivitySetschenow.hpp>
#include <Reaktor/Activity/GaseousActivity.hpp>
#include <Reaktor/Activity/GaseousActivityDuanMollerWeare.hpp>
#include <Reaktor/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktor/Activity/GaseousActivityIdeal.hpp>
#include <Reaktor/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherReed.hpp>
#include <Reaktor/Activity/MineralActivity.hpp>
#include <Reaktor/Activity/MineralActivityIdeal.hpp>

#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Functions.hpp>
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/OptimizationUtils.hpp>
#include <Reaktor/Common/Optional.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Common/ReactionEquation.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Common/StringUtils.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/TimeUtils.hpp>

#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/CoreUtils.hpp>
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Core/Reactions.hpp>
#include <Reaktor/Core/ReactionUtils.hpp>
#include <Reaktor/Core/Species.hpp>

#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktor/Equilibrium/EquilibriumUtils.hpp>

#include <Reaktor/Interfaces/Gems.hpp>
#include <Reaktor/Interfaces/Phreeqx.hpp>

#include <Reaktor/Kinetics/KineticOptions.hpp>
#include <Reaktor/Kinetics/KineticProblem.hpp>
#include <Reaktor/Kinetics/KineticResult.hpp>
#include <Reaktor/Kinetics/KineticUtils.hpp>

#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Math/Derivatives.hpp>
#include <Reaktor/Math/LagrangeInterpolator.hpp>
#include <Reaktor/Math/Roots.hpp>
#include <Reaktor/Math/MathUtils.hpp>

#include <Reaktor/Optimization/Filter.hpp>
#include <Reaktor/Optimization/Hessian.hpp>
#include <Reaktor/Optimization/KktSolver.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/OptimumSolver.hpp>
#include <Reaktor/Optimization/OptimumSolverIpfeasible.hpp>
#include <Reaktor/Optimization/OptimumSolverIpnewton.hpp>
#include <Reaktor/Optimization/OptimumSolverIpopt.hpp>
#include <Reaktor/Optimization/OptimumState.hpp>
#include <Reaktor/Optimization/Utils.hpp>

#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/BaseSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Species/ThermoParams.hpp>

#include <Reaktor/Thermodynamics/AqueousElectroState.hpp>
#include <Reaktor/Thermodynamics/AqueousElectroStateHKF.hpp>
#include <Reaktor/Thermodynamics/ThermoState.hpp>
#include <Reaktor/Thermodynamics/ThermoStateHKF.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>
#include <Reaktor/Thermodynamics/WaterElectroState.hpp>
#include <Reaktor/Thermodynamics/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzState.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzStateHGK.hpp>
#include <Reaktor/Thermodynamics/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktor/Thermodynamics/WaterThermoState.hpp>
#include <Reaktor/Thermodynamics/WaterThermoStateUtils.hpp>
#include <Reaktor/Thermodynamics/WaterUtils.hpp>
