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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SurfaceAreaModel.hpp>

namespace Reaktoro {

// Forward declarations
class PhaseList;
class Surface;

/// The function type for the generation of surfaces with given phases in the chemical system.
/// @param phases The phases composing the chemical system.
using SurfaceGenerator = Fn<Vec<Surface>(PhaseList const& phases)>;

/// Used to define a general surface.
/// @ingroup Core
class GeneralSurface
{
public:
    /// Construct a default GeneralSurface object.
    GeneralSurface();

    /// Construct a GeneralSurface object with given unique surface name.
    explicit GeneralSurface(String const& name);

    /// Set the unique name of the surface.
    auto setName(String const& name) -> GeneralSurface&;

    /// Set the area model of the surface.
    auto setAreaModel(SurfaceAreaModel const& model) -> GeneralSurface&;

    /// Set the area model of the surface (equivalent to GeneralSurface::setAreaModel).
    auto set(SurfaceAreaModel const& model) -> GeneralSurface&;

    /// Return the name of the surface.
    auto name() const -> String const&;

    /// Return the area model of the surface.
    auto areaModel() const -> SurfaceAreaModel const&;

    /// Convert this GeneralSurface object into a Surface object.
    auto operator()(PhaseList const& phases) const -> Surface;

private:
    /// The name of the surface.
    String surface_name;

    /// The area model of the surface.
    SurfaceAreaModel area_model;
};

/// Used to represent a collection of surfaces across which chemical reactions take place.
/// @ingroup Core
class Surfaces
{
public:
    /// Construct a Surfaces object.
    Surfaces();

    /// Construct a Surfaces object with given Surface, GeneralSurface, or SurfaceGenerator objects.
    /// @param surfaces The objects of type Surface, GeneralSurface, or SurfaceGenerator.
    template<typename... SurfaceConvertible>
    explicit Surfaces(SurfaceConvertible const&... surfaces)
    {
        static_assert(sizeof...(surfaces) > 0);
        addAux(surfaces...);
    }

    /// Add a surface generator into the Surfaces container.
    template<typename T>
    auto add(T const& item) -> void
    {
        static_assert(
            isConvertible<T, Surface> ||
            isConvertible<T, GeneralSurface> ||
            isConvertible<T, SurfaceGenerator>);

        if constexpr(isConvertible<T, Surface>)
        {
            surface_generators.push_back([=](PhaseList const& phases) -> Vec<Surface> { return { item }; }); // item is a Surface object
        }
        else if constexpr(isConvertible<T, GeneralSurface>)
        {
            surface_generators.push_back([=](PhaseList const& phases) -> Vec<Surface> { return { item(phases) }; }); // item is a GeneralSurface object; use operator()(PhaseList) to convert to Surface
        }
        else
        {
            surface_generators.push_back(item); // item is already a SurfaceGenerator
        }
    }

    /// Convert this Surfaces object into a vector of Surface objects.
    /// @param phases The phases composing the chemical system.
    auto convert(PhaseList const& phases) const -> Vec<Surface>;

private:
    /// The SurfaceGenerator objects collected so far with each call to Surfaces::add method.
    Vec<SurfaceGenerator> surface_generators;

    /// Add one or more SurfaceGenerator or Surface objects into the Surfaces container.
    template<typename Arg, typename... Args>
    auto addAux(const Arg& arg, const Args&... args) -> void
    {
        add(arg);
        if constexpr (sizeof...(Args) > 0)
            addAux(args...);
    }
};

} // namespace Reaktoro
