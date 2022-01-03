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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Extensions/Nasa/NasaSpecies.hpp>

namespace Reaktoro {

/// Used to parse thermodynamic database files in NASA official format.
class NasaDatabaseParser
{
public:
    /// Construct a default NasaDatabaseParser object.
    NasaDatabaseParser();

    /// Construct a copy of a NasaDatabaseParser instance.
    NasaDatabaseParser(const NasaDatabaseParser& other);

    /// Destroy this NasaDatabaseParser instance.
    virtual ~NasaDatabaseParser();

    /// Assign a NasaDatabaseParser instance to this.
    auto operator=(NasaDatabaseParser other) -> NasaDatabaseParser&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
