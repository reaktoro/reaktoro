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

#include "EquilibriumSensitivity.hpp"

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>

namespace Reaktoro {

EquilibriumSensitivity::EquilibriumSensitivity()
{}

EquilibriumSensitivity::EquilibriumSensitivity(EquilibriumSpecs const& specs)
{
    initialize(specs);
}

auto EquilibriumSensitivity::initialize(EquilibriumSpecs const& specs) -> void
{
    msystem = specs.system();
    minputs = specs.inputs();

    const EquilibriumDims dims(specs);

    const auto Nw = dims.Nw;
    const auto Nn = dims.Nn;
    const auto Np = dims.Np;
    const auto Nq = dims.Nq;
    const auto Nc = dims.Nc;

    mdndw.resize(Nn, Nw);
    mdpdw.resize(Np, Nw);
    mdqdw.resize(Nq, Nw);
    mdndc.resize(Nn, Nc);
    mdpdc.resize(Np, Nc);
    mdqdc.resize(Nq, Nc);
}

auto EquilibriumSensitivity::dndw(String const& wid) const -> VectorXdConstRef
{
    const auto idx = index(minputs, wid);
    return mdndw.col(idx);
}

auto EquilibriumSensitivity::dndw() const -> MatrixXdConstRef
{
    return mdndw;
}

auto EquilibriumSensitivity::dndw(MatrixXdConstRef data) -> void
{
    errorif(mdndw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dndw(MatrixXdConstRef).");
    errorif(mdndw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dndw(MatrixXdConstRef).");
    mdndw = data;
}

auto EquilibriumSensitivity::dpdw(String const& wid) const -> VectorXdConstRef
{
    const auto idx = index(minputs, wid);
    return mdpdw.col(idx);
}

auto EquilibriumSensitivity::dpdw() const -> MatrixXdConstRef
{
    return mdpdw;
}

auto EquilibriumSensitivity::dpdw(MatrixXdConstRef data) -> void
{
    errorif(mdpdw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dpdw(MatrixXdConstRef).");
    errorif(mdpdw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dpdw(MatrixXdConstRef).");
    mdpdw = data;
}

auto EquilibriumSensitivity::dqdw(String const& wid) const -> VectorXdConstRef
{
    const auto idx = index(minputs, wid);
    return mdqdw.col(idx);
}

auto EquilibriumSensitivity::dqdw() const -> MatrixXdConstRef
{
    return mdqdw;
}

auto EquilibriumSensitivity::dqdw(MatrixXdConstRef data) -> void
{
    errorif(mdqdw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dqdw(MatrixXdConstRef).");
    errorif(mdqdw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dqdw(MatrixXdConstRef).");
    mdqdw = data;
}

auto EquilibriumSensitivity::dndc() const -> MatrixXdConstRef
{
    return mdndc;
}

auto EquilibriumSensitivity::dndc(MatrixXdConstRef data) -> void
{
    errorif(mdndc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dndc(MatrixXdConstRef).");
    errorif(mdndc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dndc(MatrixXdConstRef).");
    mdndc = data;
}

auto EquilibriumSensitivity::dpdc() const -> MatrixXdConstRef
{
    return mdpdc;
}

auto EquilibriumSensitivity::dpdc(MatrixXdConstRef data) -> void
{
    errorif(mdpdc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dpdc(MatrixXdConstRef).");
    errorif(mdpdc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dpdc(MatrixXdConstRef).");
    mdpdc = data;
}

auto EquilibriumSensitivity::dqdc() const -> MatrixXdConstRef
{
    return mdqdc;
}

auto EquilibriumSensitivity::dqdc(MatrixXdConstRef data) -> void
{
    errorif(mdqdc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dqdc(MatrixXdConstRef).");
    errorif(mdqdc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dqdc(MatrixXdConstRef).");
    mdqdc = data;
}

auto EquilibriumSensitivity::dudw() const -> MatrixXdConstRef
{
    return mdudw;
}

auto EquilibriumSensitivity::dudc() const -> MatrixXdConstRef
{
    return mdudc;
}

auto EquilibriumSensitivity::dudw(MatrixXdConstRef data) -> void
{
    mdudw = data;
}

auto EquilibriumSensitivity::dudc(MatrixXdConstRef data) -> void
{
    mdudc = data;
}

} // namespace Reaktoro
