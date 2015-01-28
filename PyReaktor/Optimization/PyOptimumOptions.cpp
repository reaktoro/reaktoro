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

#include "PyOptimumOptions.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>

namespace Reaktor {

auto export_OptimumOptions() -> void
{
    py::class_<OptimumParamsIpnewton>("OptimumParamsIpnewton")
        .def_readwrite("mu", &OptimumParamsIpnewton::mu)
        .def_readwrite("tau", &OptimumParamsIpnewton::tau)
        .def_readwrite("mux", &OptimumParamsIpnewton::mux)
        .def_readwrite("scaling", &OptimumParamsIpnewton::scaling)
        .def_readwrite("uniform_newton_step", &OptimumParamsIpnewton::uniform_newton_step)
        ;

    py::class_<OptimumParamsIpopt>("OptimumParamsIpopt")
        .def_readwrite("mu", &OptimumParamsIpopt::mu)
        .def_readwrite("delta", &OptimumParamsIpopt::delta)
        .def_readwrite("eta_phi", &OptimumParamsIpopt::eta_phi)
        .def_readwrite("gamma_alpha", &OptimumParamsIpopt::gamma_alpha)
        .def_readwrite("gamma_phi", &OptimumParamsIpopt::gamma_phi)
        .def_readwrite("gamma_theta", &OptimumParamsIpopt::gamma_theta)
        .def_readwrite("kappa_epsilon", &OptimumParamsIpopt::kappa_epsilon)
        .def_readwrite("kappa_mu", &OptimumParamsIpopt::kappa_mu)
        .def_readwrite("kappa_sigma", &OptimumParamsIpopt::kappa_sigma)
        .def_readwrite("kappa_soc", &OptimumParamsIpopt::kappa_soc)
        .def_readwrite("s_phi", &OptimumParamsIpopt::s_phi)
        .def_readwrite("s_theta", &OptimumParamsIpopt::s_theta)
        .def_readwrite("tau_min", &OptimumParamsIpopt::tau_min)
        .def_readwrite("theta_mu", &OptimumParamsIpopt::theta_mu)
        .def_readwrite("max_iters_soc", &OptimumParamsIpopt::max_iters_soc)
        .def_readwrite("soc", &OptimumParamsIpopt::soc)
        .def_readwrite("mux", &OptimumParamsIpopt::mux)
        .def_readwrite("scaling", &OptimumParamsIpopt::scaling)
        ;

    py::enum_<OptimumMethod>("OptimumMethod")
        .value("Ipnewton", OptimumMethod::Ipnewton)
        .value("Ipopt", OptimumMethod::Ipopt)
        ;

    py::class_<OptimumOptions>("OptimumOptions")
        .def_readwrite("tolerance", &OptimumOptions::tolerance)
        .def_readwrite("max_iterations", &OptimumOptions::max_iterations)
        .def_readwrite("method", &OptimumOptions::method)
        .def_readwrite("output", &OptimumOptions::output)
        .def_readwrite("ipopt", &OptimumOptions::ipopt)
        .def_readwrite("ipnewton", &OptimumOptions::ipnewton)
        ;
}

} // namespace Reaktor
