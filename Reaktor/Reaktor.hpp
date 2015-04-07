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

#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Common/ElementUtils.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Functions.hpp>
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/InterpolationUtils.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/OptimizationUtils.hpp>
#include <Reaktor/Common/Optional.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Common/ParseUtils.hpp>
#include <Reaktor/Common/ReactionEquation.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Common/StringUtils.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/ChemicalModels.hpp>
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Core/ReactionSystem.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Utils.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktor/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktor/Interfaces/Gems.hpp>
#include <Reaktor/Interfaces/internal/PhreeqcUtils.hpp>
#include <Reaktor/Interfaces/Phreeqx.hpp>
#include <Reaktor/Kinetics/KineticOptions.hpp>
#include <Reaktor/Kinetics/KineticProblem.hpp>
#include <Reaktor/Kinetics/KineticResult.hpp>
#include <Reaktor/Kinetics/KineticSolver.hpp>
#include <Reaktor/Kinetics/KineticUtils.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Math/Derivatives.hpp>
#include <Reaktor/Math/LagrangeInterpolator.hpp>
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Math/ODE.hpp>
#include <Reaktor/Math/Roots.hpp>
#include <Reaktor/Optimization/Filter.hpp>
#include <Reaktor/Optimization/Hessian.hpp>
#include <Reaktor/Optimization/Jacobian.hpp>
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
#include <Reaktor/Thermodynamics/Activity/AqueousActivity.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityDrummond.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityDuanSun.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityHKF.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityIdeal.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityPitzer.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivityRumpf.hpp>
#include <Reaktor/Thermodynamics/Activity/AqueousActivitySetschenow.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivity.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityDuanMollerWeare.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityIdeal.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivitySpycherReed.hpp>
#include <Reaktor/Thermodynamics/Activity/MineralActivity.hpp>
#include <Reaktor/Thermodynamics/Activity/MineralActivityIdeal.hpp>
#include <Reaktor/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktor/Thermodynamics/Core/Database.hpp>
#include <Reaktor/Thermodynamics/Core/Thermo.hpp>
#include <Reaktor/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktor/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktor/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include <Reaktor/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesElectroState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>
#include <Reaktor/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktor/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktor/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktor/Thermodynamics/Reactions/MineralCatalyst.hpp>
#include <Reaktor/Thermodynamics/Reactions/MineralMechanism.hpp>
#include <Reaktor/Thermodynamics/Reactions/MineralReaction.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/ThermoData.hpp>
#include <Reaktor/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktor/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktor/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktor/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktor/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktor/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktor/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktor/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktor/Thermodynamics/Water/WaterUtils.hpp>
#include <Reaktor/Utils/GeoUtils.hpp>
#include <Reaktor/Utils/Rock.hpp>
