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

#pragma once

#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModel.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDrummondCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelDuanSunCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelRumpfCO2.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModelSetschenow.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Thermo.hpp>
#include <Reaktoro/Thermodynamics/EOS/CubicEOS.hpp>
#include <Reaktoro/Extensions/Geochemistry/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include <Reaktoro/Models/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Models/ActivityModelHKF.hpp>
#include <Reaktoro/Models/ActivityModelIdeal.hpp>
#include <Reaktoro/Models/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Models/ActivityModelCubicEOS.hpp>
#include <Reaktoro/Models/ActivityModelIdeal.hpp>
#include <Reaktoro/Models/ActivityModelSpycherPruessEnnis.hpp>
#include <Reaktoro/Models/ActivityModelSpycherReed.hpp>
#include <Reaktoro/Models/ActivityModelIdeal.hpp>
#include <Reaktoro/Models/ActivityModelRedlichKister.hpp>
#include <Reaktoro/Models/PhaseChemicalModel.hpp>
#include <Reaktoro/Models/PhaseThermoModel.hpp>
#include <Reaktoro/Models/SpeciesElectroState.hpp>
#include <Reaktoro/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktoro/Models/SpeciesThermoState.hpp>
#include <Reaktoro/Models/SpeciesThermoStateHKF.hpp>
#include <Reaktoro/Models/SpeciesElectroState.hpp>
#include <Reaktoro/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktoro/Models/SpeciesThermoState.hpp>
#include <Reaktoro/Models/SpeciesThermoStateHKF.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/LiquidPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Extensions/Geochemistry/MineralCatalyst.hpp>
#include <Reaktoro/Extensions/Geochemistry/MineralMechanism.hpp>
#include <Reaktoro/Extensions/Geochemistry/MineralReaction.hpp>
#include <Reaktoro/Extensions/Water/WaterConstants.hpp>
#include <Reaktoro/Extensions/Water/WaterElectroState.hpp>
#include <Reaktoro/Extensions/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Extensions/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Extensions/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Extensions/Water/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktoro/Extensions/Water/WaterThermoState.hpp>
#include <Reaktoro/Extensions/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Extensions/Water/WaterUtils.hpp>
