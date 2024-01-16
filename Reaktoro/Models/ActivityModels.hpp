// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Models/ActivityModels/Support/CubicEOS.hpp>
#include <Reaktoro/Models/ActivityModels/Support/IonExchangeSurface.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelCubicEOS.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDavies.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDrummond.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDuanSun.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelExtendedUNIQUAC.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelHKF.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealAqueous.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealGas.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealIonExchange.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIonExchange.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPengRobinsonPhreeqcOriginal.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPhreeqc.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPhreeqcIonicStrengthPressureCorrection.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzer.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelRedlichKister.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelRumpf.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelSetschenow.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelSpycherPruessEnnis.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelSpycherReed.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelVanLaar.hpp>

/// @defgroup ActivityModels Activity Models
/// @ingroup Models
/// @brief The module in Reaktoro in which activity models for phases are implemented.
