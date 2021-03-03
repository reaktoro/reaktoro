// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>

namespace Reaktoro {
namespace {

/// Create the lambda function that extracts temperature from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// temperature in unknown, then temperature should be extracted from vector `p`.
/// Otherwise, temperature is known and available in the Params object `w`.
/// @param specs The specifications of the equilibrium solver
auto createTemperatureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, const Params&)>
{
    const auto unknownT = specs.isTemperatureUnknown();

    if(unknownT)
        return [](VectorXrConstRef p, const Params& w) { return p[0]; };

    return [](VectorXrConstRef p, const Params& w) { return w.get("T").value(); };
}

/// Create the lambda function that extracts pressure from either `p` or `w`.
/// When the specifications of the equilibrium solver (`specs`) indicates that
/// pressure in unknown, then pressure should be extracted from vector `p`.
/// Otherwise, pressure is known and available in the Params object `w`.
/// @param specs The specifications of the equilibrium solver
auto createPressureGetterFn(const EquilibriumSpecs& specs) -> Fn<real(VectorXrConstRef, const Params&)>
{
    const auto unknownT = specs.isTemperatureUnknown();
    const auto unknownP = specs.isPressureUnknown();

    if(unknownT && unknownP)
        return [](VectorXrConstRef p, const Params& params) { return p[1]; };

    if(unknownP)
        return [](VectorXrConstRef p, const Params& params) { return p[0]; };

    return [](VectorXrConstRef p, const Params& params) { return params.get("P").value(); };
}

/// Determine the index of the variable in `(n, p, w)` that has been seeded by autodiff.
auto indexOfSeededVariable(VectorXrConstRef n, VectorXrConstRef p, const Params& w)
{
    using autodiff::grad;

    // Check for a variable v in n so that grad(v) == 1 (v has been seeded!)
    auto offset = 0;
    for(auto j = 0; j < n.size(); ++j)
        if(grad(n[j]) == 1.0)
            return j;

    // Check for a variable v in p so that grad(v) == 1 (v has been seeded!)
    offset += n.size();
    for(auto j = 0; j < p.size(); ++j)
        if(grad(p[j]) == 1.0)
            return offset + j;

    // Check for a variable v in w so that grad(v) == 1 (v has been seeded!)
    offset += p.size();
    for(auto j = 0; j < w.size(); ++j)
        if(grad(w[j].value()) == 1.0)
            return offset + j;

    // There was no seeded variable in (n, p, w)
    offset += w.size();
    return offset;
}

} // namespace

struct EquilibriumProps::Impl
{
    /// The chemical properties of the system computed at given *(n, p, w)*.
    ChemicalProps props;

    /// The dimension variables in a chemical equilibrium problem specification.
    const EquilibriumDims dims;

    /// The temperature getter function for the given equilibrium specifications. @see createTemperatureGetterFn
    const Fn<real(VectorXrConstRef, const Params&)> getT;

    /// The pressure getter function for the given equilibrium specifications. @see createPressureGetterFn
    const Fn<real(VectorXrConstRef, const Params&)> getP;

    /// The partial derivatives of the serialized chemical properties *u* with respect to *(n, p, w)*.
    MatrixXd dudnpw;

    /// The array stream used during serialize and deserialize of chemical properties.
    ArrayStream<real> stream;

    /// The flag indicating if the full Jacobian matrix is been constructed.
    bool assemblying_jacobian = false;

    /// Construct an EquilibriumProps::Impl object.
    Impl(const EquilibriumSpecs& specs)
    : props(specs.system()), dims(specs),
      getT(createTemperatureGetterFn(specs)), getP(createPressureGetterFn(specs))
    {}

    /// Update the chemical properties of the chemical system.
    auto update(VectorXrConstRef n, VectorXrConstRef p, const Params& w) -> void
    {
        const auto T = getT(p, w);
        const auto P = getP(p, w);

        props.update(T, P, n);

        if(assemblying_jacobian)
        {
            props.serialize(stream);
            const auto Nnpw = dims.Nn + dims.Np + dims.Nw;
            const auto Nu = stream.data().rows();
            const auto idx = indexOfSeededVariable(n, p, w);
            errorif(idx >= Nnpw, "Expecting an autodiff seeded variable in (n, p, w) "
                "when assemblying the full Jacobian matrix of the chemical properties.");
            dudnpw.resize(Nu, Nnpw);
            auto col = dudnpw.col(idx);
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

auto EquilibriumProps::update(VectorXrConstRef n, VectorXrConstRef p, const Params& w) -> void
{
    pimpl->update(n, p, w);
}

auto EquilibriumProps::assembleFullJacobianBegin() -> void
{
    pimpl->assemblying_jacobian = true;
}

auto EquilibriumProps::assembleFullJacobianEnd() -> void
{
    pimpl->assemblying_jacobian = false;
}

auto EquilibriumProps::chemicalProps() const -> const ChemicalProps&
{
    return pimpl->props;
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
