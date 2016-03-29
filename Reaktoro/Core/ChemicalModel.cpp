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

#include "ChemicalModel.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

ChemicalModelResult::ChemicalModelResult()
{}

ChemicalModelResult::ChemicalModelResult(unsigned nspecies)
: num_species(nspecies)
{}

auto ChemicalModelResult::operator=(const std::vector<PhaseChemicalModelResult>& results_) -> ChemicalModelResult&
{
    results = results_;
    return *this;
}

auto ChemicalModelResult::lnActivityConstants() const -> ThermoVector
{
    ThermoVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, size) = results[i].ln_activity_constants;
        offset += size;
    }
    return res;
}

auto ChemicalModelResult::lnActivityCoefficients() const -> ChemicalVector
{
    ChemicalVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, offset, size, size) = results[i].ln_activity_coefficients;
        offset += size;
    }
    return res;
}

auto ChemicalModelResult::lnActivities() const -> ChemicalVector
{
    ChemicalVector res(num_species);
    unsigned offset = 0;
    for(unsigned i = 0; i < results.size(); ++i)
    {
        const unsigned size = results[i].num_species;
        res.rows(offset, offset, size, size) = results[i].ln_activities;
        offset += size;
    }
    return res;
}

} // namespace Reaktoro
