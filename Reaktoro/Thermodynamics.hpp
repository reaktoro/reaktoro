// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelRumpf.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDrummond.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelSetschenow.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDuanSun.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealGas.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealAqueous.hpp>
#include <Reaktoro/Thermodynamics/Solids/ActivityModelVanLaar.hpp>
#include <Reaktoro/Thermodynamics/Solids/ActivityModelRedlichKister.hpp>
#include <Reaktoro/Thermodynamics/Fluids/CubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Fluids/PhaseIdentification.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelCubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelSpycherPruessEnnis.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelSpycherReed.hpp>
