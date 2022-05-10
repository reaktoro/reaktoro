// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

const auto return_internal_ref = py::return_value_policy::reference_internal;

PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXr);
PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXrRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXrConstRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXr);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXrRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXrConstRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::Indices);
