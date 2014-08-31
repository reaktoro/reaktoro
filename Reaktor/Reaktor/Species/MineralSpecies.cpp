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

#include "MineralSpecies.hpp"

// Reaktor includes
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Thermodynamics/ThermoUtils.hpp>

namespace Reaktor {

MineralSpecies::MineralSpecies()
{}

auto MineralSpecies::setMolarVolume(units::MolarVolume value) -> void
{
    molar_volume$ = value;
}

auto MineralSpecies::setThermoData(const MineralThermoData& thermoData) -> void
{
    thermo_data$ = thermoData;
}

auto MineralSpecies::density() const -> units::Density
{
    return molarMass()/molarVolume();
}

auto MineralSpecies::molarVolume() const -> units::MolarVolume
{
    return molar_volume$;
}

auto MineralSpecies::thermoData() const -> const MineralThermoData&
{
    return thermo_data$;
}

auto operator<<(std::ostream& out, const MineralSpecies& species) -> std::ostream&
{
    // Get the HKF thermodynamic data of the species
    const MineralThermoDataHKF& hkf = species.thermoData().hkf.get();

    out << static_cast<GeneralSpecies>(species);
    out << "  molar volume: " << species.molarVolume() << " m3/mol" << std::endl;
    out << "  thermo data (HKF)" << std::endl;
    out << "    Gf: " << hkf.Gf << std::endl;
    out << "    Hf: " << hkf.Hf << std::endl;
    out << "    Sr: " << hkf.Sr << std::endl;
    out << "    Vr: " << hkf.Vr << std::endl;

    if(hkf.nptrans == 0)
    {
        out << "    a: " << hkf.a.front() << std::endl;
        out << "    b: " << hkf.b.front() << std::endl;
        out << "    c: " << hkf.c.front() << std::endl;
        out << "    Tmax: " << hkf.Tmax << " K" << std::endl;
    }
    else for(int i = 0; i <= hkf.nptrans; ++i)
    {
        out << "    a[" << i << "]: " << hkf.a[i] << std::endl;
        out << "    b[" << i << "]: " << hkf.b[i] << std::endl;
        out << "    c[" << i << "]: " << hkf.c[i] << std::endl;

        if(i < hkf.nptrans)
        {
            out << "    Ttr[" << i << "]: " << hkf.Ttr[i] << std::endl;
            out << "    Htr[" << i << "]: " << hkf.Htr[i] << std::endl;
            out << "    Vtr[" << i << "]: " << hkf.Vtr[i] << std::endl;
            out << "    dPdTtr[" << i << "]: " << hkf.dPdTtr[i] << std::endl;
        }
    }

    return out;
}

} // namespace Reaktor

