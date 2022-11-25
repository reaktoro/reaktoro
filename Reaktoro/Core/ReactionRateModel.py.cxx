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
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Model.py.hxx>
#include <Reaktoro/Core/ReactionRateModel.hpp>
using namespace Reaktoro;

// Return the type names of the arguments in a Python callable object.
auto getFunctionArgumentTypes(py::object obj) -> Strings
{
    auto inspect = py::module::import("inspect");
    auto signature = inspect.attr("signature");
    auto s = signature(obj); // this will fail if obj is not a callable object
    auto params = s.attr("parameters");
    Strings res;
    auto values = py::list(params.attr("values"));
    for(auto val : values)
        res.push_back(val.attr("annotation").attr("__name__").cast<String>());
    return res;
}

// Return the type name of the first argument in a Python callable object.
// The object must support only a single argument.
auto getFunctionFirstAndOnlyArgumentType(py::object obj) -> String
{
    auto inspect = py::module::import("inspect");
    auto signature = inspect.attr("signature");
    auto s = signature(obj); // this will fail if obj is not a callable object
    auto params = s.attr("parameters");
    auto numparams = py::len(params);
    errorif(numparams != 1, "Function expected to have only one argument.");
    auto values = py::list(params.attr("values")());
    auto argtype = values[0].attr("annotation").attr("__name__").cast<String>();
    return argtype;
}

// Check if a Python object is of type ReactionRateModel.
auto isReactionRateModel(py::object obj) -> bool
{
    const auto type = obj.get_type().attr("__name__").cast<String>();
    return type == "ReactionRateModel";
}

// Check if a Python object is a callable that can be converted to a ReactionRateModel.
// If the Python object is a callable with a single argument of explicit type
// ChemicalProps, then we assume the Python object can be safely converted to a
// ReactionRateModel. If this is not the case, a runtime error will happen at
// the time the converted model is evaluated for the first time.
auto isReactionRateModelConvertible(py::object obj) -> bool
{
    // Check if obj is a Python callable (function or function object with
    // __call__ method) with a single argument of explicit type ChemicalProps
    try { return getFunctionFirstAndOnlyArgumentType(obj) == "ChemicalProps"; }
    catch(...) { return false; }
}

// Check if a Python object is of type ReactionRateModelGenerator.
auto isReactionRateModelGenerator(py::object obj) -> bool
{
    // Determining whether a Python object is a ReactionRateModelGenerator
    // object, created by Reaktoro from C++, we need to check if the PyCapsule
    // wrapper created by pybind11 has any hint about
    // ReactionRateModelGenerator. After checking all attributes of such a
    // PyCapsule object, it was found that the __doc__ attribute has the
    // information we want. When this is the case, obj.__doc__ contains a string
    // such as this:
    //
    // '(arg0: reaktoro.reaktoro4py.ReactionRateModelGeneratorArgs) -> reaktoro.reaktoro4py.ReactionRateModel\n'
    //
    // we use this then (no better alternative known at the moment!) to
    // determine if obj is a ReactionRateModelGenerator object.
    try {
        const auto doc = obj.attr("__doc__").cast<String>();
        return doc.find("ReactionRateModelGeneratorArgs") != String::npos;
    } catch(...) {
        return false; // in case obj.__doc__ is None and cannot be converted to string
    }
}

// Check if a Python object is a callable that can be converted to a ReactionRateModelGenerator.
// If the Python object is a callable with a single argument of explicit type
// ReactionRateModelGeneratorArgs, then we assume the Python object can be
// safely converted to a ReactionRateModelGenerator. If this is not the case, a
// runtime error will happen at the time the converted model is evaluated for
// the first time.
auto isReactionRateModelGeneratorConvertible(py::object obj) -> bool
{
    // Check if obj is a Python callable (function or function object with
    // __call__ method) with a single argument of explicit type ReactionRateModelGeneratorArgs
    try { return getFunctionFirstAndOnlyArgumentType(obj) == "ReactionRateModelGeneratorArgs"; }
    catch(...) { return false; }
}

// Create a ReactionRateModel object from a Python callable object.
// This method expects that the Python callable has only one argument of type
// ChemicalProps. This conversion function is implemented in a way that allow
// users to define reaction rate models from Python that return real instead of
// a ReactionRate object. If return type is real, this is converted to a
// ReactionRate object.
auto createReactionRateModel(py::object obj) -> ReactionRateModel
{
    errorif(!isReactionRateModelConvertible(obj), "Expecting a Python callable object with just one argument of type ChemicalProps.");

    return [=](ChemicalProps const& props) -> ReactionRate {
        auto res = obj(props);
        try { return ReactionRate(res.cast<real>()); }
        catch(...) {
            try { return res.cast<ReactionRate>(); }
            catch(...) {
                errorif(true, "Your reaction rate model definition does not return a value convertible to autodiff.real or ReactionRate.");
                return {};
            }
        }
    };
}

// Create a ReactionRateModelGenerator object from a Python callable object.
// This method expects that the Python callable has only one argument of type
// ReactionRateModelGeneratorArgs.
auto createReactionRateModelGenerator(py::object obj) -> ReactionRateModelGenerator
{
    errorif(!isReactionRateModelGeneratorConvertible(obj), "Expecting a Python callable object with just one argument of type ReactionRateModelGeneratorArgs.");

    return [=](ReactionRateModelGeneratorArgs args) -> ReactionRateModel
    {
        auto resfn = obj(args);
        return createReactionRateModel(resfn);
    };
}

void exportReactionRateModel(py::module& m)
{
    auto cls = exportModelMethodsOnly<ReactionRate, ChemicalProps const&>(m, "ReactionRateModel")
        .def(py::init<>())
        .def(py::init(&createReactionRateModel)) // Enable construction of a ReactionRateModel object using a Python callable object with single argument of type ChemicalProps
        ;

    py::implicitly_convertible<Fn<real(ChemicalProps const&)>, ReactionRateModel>();
    py::implicitly_convertible<Fn<ReactionRate(ChemicalProps const&)>, ReactionRateModel>();

    py::class_<ReactionRateModelGeneratorArgs>(m, "ReactionRateModelGeneratorArgs")
        .def(py::init<String const&, ReactionEquation const&, Database const&, SpeciesList const&, PhaseList const&, SurfaceList const&>())
        .def_property_readonly("name", [](ReactionRateModelGeneratorArgs const& self) { return self.name; }, "The name of the reaction for which the rate model is generated.")
        .def_property_readonly("equation", [](ReactionRateModelGeneratorArgs const& self) { return self.equation; }, "The equation of the reaction for which the rate model is generated.")
        .def_property_readonly("database", [](ReactionRateModelGeneratorArgs const& self) { return self.database; }, "The thermodynamic database used to construct the chemical system where the reaction belongs to.")
        .def_property_readonly("species", [](ReactionRateModelGeneratorArgs const& self) { return self.species; }, "The species in the chemical system where the reaction belongs to.")
        .def_property_readonly("phases", [](ReactionRateModelGeneratorArgs const& self) { return self.phases; }, "The phases in the chemical system where the reaction belongs to.")
        .def_property_readonly("surfaces", [](ReactionRateModelGeneratorArgs const& self) { return self.surfaces; }, "The surfaces in the chemical system where the reaction belongs to.")
        ;
}
