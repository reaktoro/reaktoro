// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumOptions.hpp>

namespace Reaktoro {

auto export_OptimumOptions() -> void
{
    py::class_<OptimumParamsIpNewton>("OptimumParamsIpNewton")
        .def_readwrite("mu", &OptimumParamsIpNewton::mu)
        .def_readwrite("tau", &OptimumParamsIpNewton::tau)
        .def_readwrite("mux", &OptimumParamsIpNewton::mux)
        .def_readwrite("scaling", &OptimumParamsIpNewton::scaling)
        .def_readwrite("uniform_newton_step", &OptimumParamsIpNewton::uniform_newton_step)
        ;

    py::class_<OptimumParamsIpOpt>("OptimumParamsIpOpt")
        .def_readwrite("mu", &OptimumParamsIpOpt::mu)
        .def_readwrite("delta", &OptimumParamsIpOpt::delta)
        .def_readwrite("eta_phi", &OptimumParamsIpOpt::eta_phi)
        .def_readwrite("gamma_alpha", &OptimumParamsIpOpt::gamma_alpha)
        .def_readwrite("gamma_phi", &OptimumParamsIpOpt::gamma_phi)
        .def_readwrite("gamma_theta", &OptimumParamsIpOpt::gamma_theta)
        .def_readwrite("kappa_epsilon", &OptimumParamsIpOpt::kappa_epsilon)
        .def_readwrite("kappa_mu", &OptimumParamsIpOpt::kappa_mu)
        .def_readwrite("kappa_sigma", &OptimumParamsIpOpt::kappa_sigma)
        .def_readwrite("kappa_soc", &OptimumParamsIpOpt::kappa_soc)
        .def_readwrite("s_phi", &OptimumParamsIpOpt::s_phi)
        .def_readwrite("s_theta", &OptimumParamsIpOpt::s_theta)
        .def_readwrite("tau_min", &OptimumParamsIpOpt::tau_min)
        .def_readwrite("theta_mu", &OptimumParamsIpOpt::theta_mu)
        .def_readwrite("max_iters_soc", &OptimumParamsIpOpt::max_iters_soc)
        .def_readwrite("soc", &OptimumParamsIpOpt::soc)
        .def_readwrite("mux", &OptimumParamsIpOpt::mux)
        .def_readwrite("scaling", &OptimumParamsIpOpt::scaling)
        ;

    py::enum_<OptimumMethod>("OptimumMethod")
        .value("IpNewton", OptimumMethod::IpNewton)
        .value("IpOpt", OptimumMethod::IpOpt)
        ;

    py::class_<OptimumOutput, py::bases<OutputterOptions>>("OptimumOutput")
        .def_readwrite("xprefix", &OptimumOutput::xprefix)
        .def_readwrite("yprefix", &OptimumOutput::yprefix)
        .def_readwrite("zprefix", &OptimumOutput::zprefix)
        .def_readwrite("xnames", &OptimumOutput::xnames)
        .def_readwrite("ynames", &OptimumOutput::ynames)
        .def_readwrite("znames", &OptimumOutput::znames)
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

} // namespace Reaktoro
