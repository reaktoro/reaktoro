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
#include <Reaktoro/Core/Database.hpp>

// Reaktoro includes
namespace Reaktoro {

/// Used to support thermodynamic databases in NASA CEA format.
/// @ingroup Databases
class NasaDatabase : public Database
{
public:
    /// Construct a default NasaDatabase instance.
    NasaDatabase();

    /// Construct a copy of a NasaDatabase instance.
    NasaDatabase(const NasaDatabase& other);

    /// Construct a NasaDatabase instance with given species.
    explicit NasaDatabase(const Vec<Species>& species);

    /// Destroy this NasaDatabase instance.
    ~NasaDatabase();

    /// Assign another NasaDatabase instance to this.
    auto operator=(NasaDatabase other) -> NasaDatabase&;


private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};


} // namespace Reaktoro
