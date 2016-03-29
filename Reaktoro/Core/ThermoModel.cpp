// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "ThermoModel.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseThermoModel.hpp>

namespace Reaktoro {

ThermoModelResult::ThermoModelResult()
{}

ThermoModelResult::ThermoModelResult(unsigned nspecies)
: num_species(nspecies)
{}

auto ThermoModelResult::operator=(const std::vector<PhaseThermoModelResult>& results_) -> ThermoModelResult&
{
    results = results_;
    return *this;
}

auto ThermoModelResult::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].standard_partial_molar_gibbs_energies;
        offset += size;
    }
    return res;
}

auto ThermoModelResult::standardPartialMolarEnthalpies() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].standard_partial_molar_enthalpies;
        offset += size;
    }
    return res;
}

auto ThermoModelResult::standardPartialMolarVolumes() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].standard_partial_molar_volumes;
        offset += size;
    }
    return res;
}

auto ThermoModelResult::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].standard_partial_molar_heat_capacities_cp;
        offset += size;
    }
    return res;
}

auto ThermoModelResult::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].standard_partial_molar_heat_capacities_cv;
        offset += size;
    }
    return res;
}

} // namespace Reaktoro
