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

#include "Core.hpp"

// Reaktor includes
#include "Core/ChemicalSystem.hpp"
#include "Core/Element.hpp"
#include "Core/Partition.hpp"
#include "Core/Phase.hpp"
#include "Core/Reaction.hpp"
#include "Core/Species.hpp"
#include "Core/Utils.hpp"

namespace Reaktor {

auto export_Core() -> void
{
	export_Element();
	export_Species();
	export_Phase();
	export_ChemicalSystem();
}

} // namespace Reaktor
