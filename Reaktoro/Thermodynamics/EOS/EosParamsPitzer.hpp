// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <vector>
#include <string>

namespace Reaktoro {

/// A struct that defines the interaction parameters for the Pitzer model described in Harvie, Møller, and Weare (1984).
/// This struct stores the binary and ternary interaction parameters for the Pitzer model described in *Harvie, C.E., Møller,
/// N. and Weare, J.H., 1984. The prediction of mineral solubilities in natural waters: The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O
/// system to high ionic strengths at 25°C. Geochimica et Cosmochimica Acta, 48(4), pp.723–751*.
struct EosParamsPitzer
{
	/// Construct a default EosParamsPitzer instance.
	EosParamsPitzer();

	/// Construct a EosParamsPitzer instance by parsing a file with Pitzer parameters.
	EosParamsPitzer(std::string filename);

    /// The binary interaction parameters \f$ \beta_0 \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::vector<double>>> beta0;

    /// The binary interaction parameters \f$ \beta_1 \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::vector<double>>> beta1;

    /// The binary interaction parameters \f$ \beta_2 \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::vector<double>>> beta2;

    /// The binary interaction parameters \f$ \C^{\phi} \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::vector<double>>> cphi;

    /// The binary interaction parameters \f$ \theta \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, double>> theta;

    /// The binary interaction parameters \f$ \lambda \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, double>> lambda;

    /// The ternary interaction parameters \f$ \psi \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::string, double>> psi;

    /// The ternary interaction parameters \f$ \zeta \f$ of the Pizer model of Harvie, Møller, and Weare (1984).
    std::vector<std::tuple<std::string, std::string, std::string, double>> zeta;
};

} // namespace Reaktoro
