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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Transport/TransportOptions.hpp>

namespace Reaktoro {

void exportTransportOptions(py::module& m)
{
    py::enum_<FiniteVolumeMethod>(m, "FiniteVolumeMethod")
        .value("FullImplicit", FiniteVolumeMethod::FullImplicit)
        .value("ImplicitExplicit", FiniteVolumeMethod::ImplicitExplicit)
        .value("FluxLimitersImplicitExplicit", FiniteVolumeMethod::FluxLimitersImplicitExplicit)
        ;

    py::class_<TransportOptions>(m, "TransportOptions")
        .def("finite_volume_method", &TransportOptions::finite_volume_method)
        ;
}

} // namespace Reaktoro
