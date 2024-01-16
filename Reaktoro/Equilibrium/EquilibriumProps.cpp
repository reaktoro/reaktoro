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

/// Alias for a function type that extracts temperature or pressure for either variables `p` or `w`.
using PropertyGetterFn = Fn<real(VectorXrConstRef, VectorXrConstRef)>;

/// Create the lambda function that extracts temperature from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// temperature in unknown, then temperature should be extracted from vector `p`.
/// Otherwise, temperature is known and available in `w`.
/// @param specs The specifications of the equilibrium solver
auto createTemperatureGetterFn(const EquilibriumSpecs& specs) -> PropertyGetterFn
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
auto createPressureGetterFn(const EquilibriumSpecs& specs) -> PropertyGetterFn
{
    const auto iPw = index(specs.namesInputs(), "P");
    const auto iPp = index(specs.namesControlVariables(), "P");
    const auto Nw = specs.numInputs();
    if(iPw < Nw)
        return [iPw](VectorXrConstRef p, VectorXrConstRef w) { return w[iPw]; };
    return [iPp](VectorXrConstRef p, VectorXrConstRef w) { return p[iPp]; };
}

} // namespace

struct EquilibriumProps::Impl
{
    ChemicalState state;               ///< The chemical state of the system and its properties at given *(n, p, w)*.
    EquilibriumSpecs const specs;      ///< The specifications of the equilibrium problem.
    EquilibriumDims const dims;        ///< The dimension variables in a chemical equilibrium problem specification.
    PropertyGetterFn const getT;       ///< The temperature getter function for the given equilibrium specifications.
    PropertyGetterFn const getP;       ///< The pressure getter function for the given equilibrium specifications.
    MatrixXd dudnpw;                   ///< The partial derivatives of the serialized chemical properties *u* with respect to *(n, p, w)*.
    ArrayStream<real> stream;          ///< The array stream used during serialize and deserialize of chemical properties.
    bool assemblying_jacobian = false; ///< The flag indicating if the full Jacobian matrix is been constructed.

    /// Construct an EquilibriumProps::Impl object.
    Impl(const EquilibriumSpecs& specs)
    : state(specs.system()), specs(specs), dims(specs), getT(createTemperatureGetterFn(specs)), getP(createPressureGetterFn(specs))
    {
        // Initialize Jacobian matrix dudnpw with zeros (to avoid uninitialized values)
        state.props().serialize(stream);
        const auto Nu = stream.data().rows();
        const auto Nnpw = dims.Nn + dims.Np + dims.Nw;
        dudnpw = zeros(Nu, Nnpw);
    }

    /// Update the chemical properties of the chemical system.
    auto update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel) -> void
    {
        auto const T = getT(p, w);
        auto const P = getP(p, w);

        if(useIdealModel)
            state.updateIdeal(T, P, n);
        else state.update(T, P, n);
    }

    /// Update the chemical properties of the chemical system.
    auto update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel, long inpw) -> void
    {
        // Update the actual properties of the system
        update(n, p, w, useIdealModel);

        // Collect the derivatives of the chemical properties wrt some seeded variable in n, p, w.
        if(assemblying_jacobian && inpw != -1)  // inpw === -1 if seeded variable is some variable in q (the amounts of implicit titrants)
        {
            const auto Nnpw = dims.Nn + dims.Np + dims.Nw;
            assert(inpw < Nnpw);
            state.props().serialize(stream);
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

auto EquilibriumProps::update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel) -> void
{
    pimpl->update(n, p, w, useIdealModel);
}

auto EquilibriumProps::update(VectorXrConstRef n, VectorXrConstRef p, VectorXrConstRef w, bool useIdealModel, long inpw) -> void
{
    pimpl->update(n, p, w, useIdealModel, inpw);
}

auto EquilibriumProps::assembleFullJacobianBegin() -> void
{
    pimpl->assemblying_jacobian = true;
    pimpl->dudnpw.fill(0.0); // initialize with zeros to remove previous derivative values
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
