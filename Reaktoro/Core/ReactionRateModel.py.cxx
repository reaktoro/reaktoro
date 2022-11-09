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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ReactionRateModel.hpp>
#include <Reaktoro/Core/Model.py.hxx>
using namespace Reaktoro;

void exportReactionRateModel(py::module& m)
{
    // Define a function that will serve as a constructor for a ReactionRateModel object. This
    // function is special in the sesne that it allow users to define reaction rate models from
    // Python that return float instead of a ReactionRate object. This custom constructor gets this
    // float value, encapsulate it into a ReactionRate object, which is returned.
    auto createReactionRateModel = [](py::function fn)
    {
        return [=](ChemicalProps const& props) -> ReactionRate {
            auto res = fn(props);
            try { return ReactionRate(res.cast<double>()); }
            catch(...) {
                try { return res.cast<ReactionRate>(); }
                catch(...) {
                    errorif(true, "Your reaction rate definition does not return a value convertible to float or ReactionRate.");
                    return {};
                }
            }
        };
    };

    auto cls = exportModelMethodsOnly<ReactionRate, ChemicalProps const&>(m, "ReactionRateModel")
        .def(py::init<>())
        .def(py::init(createReactionRateModel))
        ;

    py::implicitly_convertible<Fn<double(ChemicalProps const&)>, ReactionRateModel>();
    py::implicitly_convertible<Fn<ReactionRate(ChemicalProps const&)>, ReactionRateModel>();
}
