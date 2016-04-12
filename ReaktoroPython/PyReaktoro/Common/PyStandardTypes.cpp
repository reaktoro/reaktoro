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

#include <PyReaktoro/Common/PyConverters.hpp>
#include "PyStandardTypes.hpp"

// PyReaktoro includes

namespace Reaktoro {

auto export_StandardTypes() -> void
{
    export_std_vector_with_str<char>("CharVector");
    export_std_vector_with_str<bool>("BoolVector");
    export_std_vector_with_str<int>("IntVector");
    export_std_vector_with_str<float>("FloatVector");
    export_std_vector_with_str<double>("DoubleVector");
    export_std_vector_with_str<std::string>("StringVector");
    export_std_vector_with_str<std::size_t>("SizetVector");
    export_std_vector<std::vector<std::size_t>>("SizetVectorVector");
    export_std_map_with_str<std::string,double>("StringDoubleMap");
}

} // namespace Reaktoro
