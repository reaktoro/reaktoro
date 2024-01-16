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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// A class that still needs to be commented.
class Class
{
public:
    /// Construct a default Class object.
    Class();

    /// Construct a copy of a Class object.
    Class(const Class& other);

    /// Destroy this Class object.
    ~Class();

    /// Assign a copy of a Class object to this.
    auto operator=(Class other) -> Class&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
