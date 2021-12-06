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

#ifndef REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP
#define REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP

// C++ includes
#include <map>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;
class EUNIQUACParams;

/// Return an equation of state for an aqueous phase based on the E-UNIQUAC activity model.
/// @param mixture The aqueous mixture instance
/// @param params The parameters for the E-UNIQUAC activity model.
/// @return The activity model function for the aqueous phase
/// @see AqueousMixture, EUNIQUACParams, PhaseChemicalModel
auto aqueousChemicalModelEUNIQUAC(const AqueousMixture& mixture, const EUNIQUACParams& params) -> PhaseChemicalModel;

class EUNIQUACParams
{
public:
    /// Construct a default EUNIQUACParams instance.
    EUNIQUACParams();

    /// Get UNIQUAC r_i parameter given a species name
    auto ri(std::string name) const -> double;

    auto ri(std::string name, double value) -> void;

    auto ri(const std::map<std::string, double>& pairs) -> void;

    /// Get UNIQUAC q_i parameter given a species name
    auto qi(std::string name) const -> double;

    auto qi(std::string name, double value) -> void;

    auto qi(const std::map<std::string, double>& pairs) -> void;

    auto setDTUvalues() -> void;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;

    /// A map to identify species indices in entalphic BIP matrices
    std::unordered_map<std::string, int> bips_species_id_map;
};

}

#endif //REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP
