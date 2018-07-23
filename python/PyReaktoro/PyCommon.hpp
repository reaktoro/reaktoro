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

// PyReaktoro includes
#include <PyReaktoro/Common/PyAutoDiff.hpp>
#include <PyReaktoro/Common/PyEigen.hpp>
#include <PyReaktoro/Common/PyMatrix.hpp>
#include <PyReaktoro/Common/PyOutputter.hpp>
#include <PyReaktoro/Common/PyReactionEquation.hpp>
#include <PyReaktoro/Common/PyStringList.hpp>
#include <PyReaktoro/Common/PyStandardTypes.hpp>
#include <PyReaktoro/Common/PyUnits.hpp>

namespace Reaktoro {

inline auto export_Common() -> void
{
    export_AutoDiff();
    export_Eigen();
    export_Matrix();
    export_ReactionEquation();
    export_StandardTypes();
    export_StringList();
    export_Outputter();
    export_Units();
}

} // namespace Reaktoro
