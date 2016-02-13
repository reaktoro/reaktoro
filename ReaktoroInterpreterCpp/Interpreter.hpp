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

#pragma once

// C++ includes
#include <istream>
#include <map>
#include <memory>
#include <string>

namespace Reaktoro {

/// A type used to define operations that interpret a Reaktoro script file.
class Interpreter
{
public:
    /// Construct a default Interpreter instance.
    Interpreter();

    /// Construct a copy of an Interpreter instance
    Interpreter(const Interpreter& other);

    /// Destroy this instance
    virtual ~Interpreter();

    /// Assign an Interpreter instance to this instance
    auto operator=(Interpreter other) -> Interpreter&;

    /// Execute a Reaktoro input script as string.
    auto execute(std::string str) -> void;

    /// Execute a Reaktoro input script as a file.
    auto execute(std::istream& stream) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
