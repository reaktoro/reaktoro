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
#include <fstream>

// json includes
#include <nlohmann/json.hpp>

namespace Reaktoro {

/// A type alias for nlohmann::json
using json = nlohmann::json;

/// Provides a convenient mechanism for outputing json files.
class JsonOutput
{
public:
    /// Construct a JsonOutput object.
    JsonOutput(std::string filename)
    : file(filename)
    {}

    /// Output a value (converted to json object) to a file.
    template<typename T>
    auto operator<<(const T& val) -> JsonOutput&
    {
        file << json(val);
        file.close();
        return *this;
    }

private:
    std::ofstream file;
};

} // namespace Reaktoro
