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

#pragma once

// Optima includes
#include <Optima/IPFilter/IPFilterOptions.hpp>

namespace Reaktor {

struct EquilibriumOptions : Optima::IPFilterOptions
{
    EquilibriumOptions();

    struct Refinement
    {
        /**
         * The maximum number of iterations of the refinement algorithm
         */
        unsigned max_iterations = 100;

        /**
         * The tolerance value used in the convergence checking of the refinement algorithm
         */
        double tolerance = 1.0e-8;
    };

    /**
     * The parameters for the refinement algorithm
     */
    Refinement refinement;
};

} // namespace Reaktor
