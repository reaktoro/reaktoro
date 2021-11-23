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

#include <Reaktoro/Thermodynamics/Water/WaterThermoProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoPropsUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzPropsWagnerPruss.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroPropsJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzPropsHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelRumpf.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDrummond.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelSetschenow.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDavies.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDuanSun.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealGas.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealAqueous.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealIonExchange.hpp>
#include <Reaktoro/Thermodynamics/Solids/ActivityModelVanLaar.hpp>
#include <Reaktoro/Thermodynamics/Solids/ActivityModelRedlichKister.hpp>
#include <Reaktoro/Thermodynamics/Fluids/CubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelCubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelSpycherPruessEnnis.hpp>
#include <Reaktoro/Thermodynamics/Fluids/ActivityModelSpycherReed.hpp>
#include <Reaktoro/Thermodynamics/Surface/ActivityModelIonExchange.hpp>
#include <Reaktoro/Thermodynamics/Surface/ActivityModelSurfaceComplexation.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeProps.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>

/// @defgroup Thermodynamics Thermodynamics
/// This is the thermodynamics module in Reaktoro, in which various thermodynamic models and related concepts are implemented.
