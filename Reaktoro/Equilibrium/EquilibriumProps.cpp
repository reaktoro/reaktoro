// Reaktoro is a unified framework for modeling chemically reactive phases.
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

#include "EquilibriumProps.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ArrayStream.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {
namespace {

/// Create the lambda function that extracts temperature from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// temperature in unknown, then temperature should be extracted from vector `p`.
/// Otherwise, temperature is known and available in `w`.
/// @param specs The specifications of the equilibrium solver
auto createTemperatureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, VectorXrConstRef)>
{
    const auto iTw = index(specs.namesInputs(), "T");
    const auto iTp = index(specs.namesControlVariables(), "T");
    const auto Nw = specs.numInputs();
    if(iTw < Nw)
        return [iTw](VectorXrConstRef p, VectorXrConstRef w) { return w[iTw]; };
    return [iTp](VectorXrConstRef p, VectorXrConstRef w) { return p[iTp]; };
}

/// Create the lambda function that extracts pressure from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// pressure in unknown, then pressure should be extracted from vector `p`.
/// Otherwise, pressure is known and available in `w`.
/// @param specs The specifications of the equilibrium solver
auto createPressureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, VectorXrConstRef)>
{
    const auto iPw = index(specs.namesInputs(), "P");
    const auto iPp = index(specs.namesControlVariables(), "P");
    const auto Nw = specs.numInputs();
    if(iPw < Nw)
        return [iPw](VectorXrConstRef p, VectorXrConstRef w) { return w[iPw]; };
    return [iPp](VectorXrConstRef p, VectorXrConstRef w) { return p[iPp]; };
}

/// Create the lambda function that extracts surface areas of reactive phase interfaces from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// a reactive surface in unknown, then the corresponding surface area should
/// be extracted from vector `p`. Otherwise, it is known and available in `w`.
/// @param specs The specifications of the equilibrium solver
auto createSurfaceAreaGetterFns(const EquilibriumSpecs& specs) -> Vec<Fn<real(VectorXrConstRef, VectorXrConstRef)>>
{
    auto const& system = specs.system();
    auto const& surfaces = system.surfaces();

    auto const num_surfaces = surfaces.size();

    auto const Nw = specs.numInputs();
    auto const Np = specs.numControlVariablesP();

    Vec<Fn<real(VectorXrConstRef, VectorXrConstRef)>> fns;

    for(auto const& surface : system.surfaces())
    {
        auto const id = "surfaceArea[" + surface.name() + "]"; // surfaceArea[AqueousPhase:GaseousPhase], surfaceArea[AqueousPhase:Quartz], surfaceArea[Calcite]
        auto const iSAw = index(specs.namesInputs(), id);
        auto const iSAp = index(specs.namesControlVariables(), id);

        errorif(iSAw >= Nw && iSAp >= Np, "Expecting surface area with name `", surface.name(), "` to be either an input or a p control variable in the equilibrium calculation.");

        if(iSAw < Nw) fns.push_back( [iSAw](VectorXrConstRef p, VectorXrConstRef w) -> real { return w[iSAw]; } );
        else fns.push_back( [iSAp](VectorXrConstRef p, VectorXrConstRef w) -> real { return p[iSAp]; } );
    }

    return fns;
}

} // namespace

struct EquilibriumProps::Impl
{
    /// The chemical state of the system and its properties at given *(n, p, w)*.
    ChemicalState state;

    /// The specifications of the equilibrium problem
    const EquilibriumSpecs specs;

    /// The dimension variables in a chemical equilibrium problem specification.
    const EquilibriumDims dims;

    /// The temperature getter function for the given equilibrium specifications. @see createTemperatureGetterFn
    const Fn<real(VectorXrConstRef, VectorXrConstRef)> getT;

    /// The pressure getter function for the given equilibrium specifications. @see createPressureGetterFn
    const Fn<real(VectorXrConstRef, VectorXrConstRef)> getP;

    /// The surface area getter functions for the given equilibrium specifications. @see createSurfaceAreaGetterFns
    const Vec<Fn<real(VectorXrConstRef, VectorXrConstRef)>> getSAs;

    /// The values of the model parameters before they are altered in the update method.
    VectorXr params0;

    /// The partial derivatives of the serialized chemical properties *u* with respect to *(n, p, w)*.
    MatrixXd dudnpw;

    /// The surface areas of the reactive phase interfaces in the system (in m2).
    VectorXr s;

    /// The array stream used during serialize and deserialize of chemical properties.
    ArrayStream<real> stream;

    /// The flag indicating if the full Jacobian matrix is been constructed.
    bool assemblying_jacobian = false;

    /// Construct an EquilibriumProps::Impl object.
    Impl(const EquilibriumSpecs& specs)
    : state(specs.system()),
      specs(specs),
      dims(specs),
      getT(createTemperatureGetterFn(specs)),
      getP(createPressureGetterFn(specs)),
      getSAs(createSurfaceAreaGetterFns(specs)),
      s(zeros(specs.system().surfaces().size()))
    {}

    /// Update the chemical properties of the chemical system.
    auto update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel, long inpw) -> void
    {
        // Get temperature and pressure of the system, either available in p or w
        const auto T = getT(p, w);
        const auto P = getP(p, w);

        // Update the surface areas of the reactive phase interfaces in the system, either available in p or w
        for(auto const& [i, getSA] : enumerate(getSAs))
            s[i] = getSA(p, w);

        // The model parameters considered inputs in the equilibrium calculation.
        auto params = specs.params();

        // The indices of the model params among the inputs.
        const auto& iparams = specs.indicesParams();

        // Store the current values of the model parameters
        params0.resize(params.size());
        for(const auto& [i, param] : enumerate(params))
            params0[i] = param.value();

        // Before updating the chemical properties, change the model parameters
        // that are input in the chemical equilibrium calculation.
        for(auto i = 0; i < params0.size(); ++i)
            params[i].value() = w[iparams[i]];

        // Perform the update of the chemical properties of the system.
        // If there were model parameters changed above, the chemical
        // properties computed below will be affected.
        if(useIdealModel)
            state.updateIdeal(T, P, n, s);
        else state.update(T, P, n, s);

        // Recover here the original state of the model parameters changed above.
        for(auto i = 0; i < params0.size(); ++i)
            params[i].value() = params0[i];

        // Collect the derivatives of the chemical properties wrt some seeded variable in n, p, w.
        if(assemblying_jacobian && inpw != -1)  // inpw === -1 if seeded variable is some variable in q (the amounts of implicit titrants)
        {
            const auto Nnpw = dims.Nn + dims.Np + dims.Nw;
            assert(inpw < Nnpw);
            state.props().serialize(stream);
            const auto Nu = stream.data().rows();
            dudnpw.resize(Nu, Nnpw);
            auto col = dudnpw.col(inpw);
            const auto size = col.size();
            for(auto i = 0; i < size; ++i)
                col[i] = grad(stream.data()[i]);
        }
    }

    /// Return the partial derivatives *du/dn*.
    auto dudn() const -> MatrixXdConstRef
    {
        return dudnpw.leftCols(dims.Nn);
    }

    /// Return the partial derivatives *du/dp*.
    auto dudp() const -> MatrixXdConstRef
    {
        return dudnpw.middleCols(dims.Nn, dims.Np);
    }

    /// Return the partial derivatives *du/dw*.
    auto dudw() const -> MatrixXdConstRef
    {
        return dudnpw.rightCols(dims.Nw);
    }
};

EquilibriumProps::EquilibriumProps(const EquilibriumSpecs& specs)
: pimpl(new Impl(specs))
{}

EquilibriumProps::EquilibriumProps(const EquilibriumProps& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumProps::~EquilibriumProps()
{}

auto EquilibriumProps::operator=(EquilibriumProps other) -> EquilibriumProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumProps::update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel, long inpw) -> void
{
    pimpl->update(n, p, w, useIdealModel, inpw);
}

auto EquilibriumProps::assembleFullJacobianBegin() -> void
{
    pimpl->assemblying_jacobian = true;
}

auto EquilibriumProps::assembleFullJacobianEnd() -> void
{
    pimpl->assemblying_jacobian = false;
}

auto EquilibriumProps::chemicalState() const -> const ChemicalState&
{
    return pimpl->state;
}

auto EquilibriumProps::chemicalProps() const -> const ChemicalProps&
{
    return pimpl->state.props();
}

auto EquilibriumProps::dudn() const -> MatrixXdConstRef
{
    return pimpl->dudn();
}

auto EquilibriumProps::dudp() const -> MatrixXdConstRef
{
    return pimpl->dudp();
}

auto EquilibriumProps::dudw() const -> MatrixXdConstRef
{
    return pimpl->dudw();
}

} // namespace Reaktoro
