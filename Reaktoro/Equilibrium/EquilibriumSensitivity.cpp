// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

struct EquilibriumSensitivity::Impl
{
    /// The chemical system associated with the sensitivity derivatives.
    ChemicalSystem system;

    /// The input parameters *c* in the chemical equilibrium problem specifications.
    Params params;

    /// The derivatives of the species amounts *n* with respect to input parameters *c*.
    MatrixXd dndc;

    /// The derivatives of the control variables *p* with respect to input parameters *c*.
    MatrixXd dpdc;

    /// The derivatives of the control variables *q* with respect to input parameters *c*.
    MatrixXd dqdc;

    /// The derivatives of the species amounts *n* with respect to component amounts *b*.
    MatrixXd dndb;

    /// The derivatives of the control variables *p* with respect to component amounts *b*.
    MatrixXd dpdb;

    /// The derivatives of the control variables *q* with respect to component amounts *b*.
    MatrixXd dqdb;

    /// Construct a default Impl object.
    Impl()
    {}

    /// Construct a default Impl object.
    Impl(const EquilibriumSpecs& specs)
    {
        initialize(specs);
    }

    /// Initialize this EquilibriumSensitivity object with given equilibrium problem specifications.
    auto initialize(const EquilibriumSpecs& specs) -> void
    {
        system = specs.system();
        params = specs.params();

        const EquilibriumDims dims(specs);

        const auto Nc = params.size();
        const auto Nn = dims.Nn;
        const auto Np = dims.Np;
        const auto Nq = dims.Nq;
        const auto Nb = dims.Nb;

        dndc.resize(Nn, Nc);
        dpdc.resize(Np, Nc);
        dqdc.resize(Nq, Nc);
        dndb.resize(Nn, Nb);
        dpdb.resize(Np, Nb);
        dqdb.resize(Nq, Nb);
    }
};

EquilibriumSensitivity::EquilibriumSensitivity()
: pimpl(new Impl())
{}

EquilibriumSensitivity::EquilibriumSensitivity(const EquilibriumSpecs& specs)
: pimpl(new Impl(specs))
{}

EquilibriumSensitivity::EquilibriumSensitivity(const EquilibriumSensitivity& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSensitivity::~EquilibriumSensitivity()
{}

auto EquilibriumSensitivity::operator=(EquilibriumSensitivity other) -> EquilibriumSensitivity&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSensitivity::initialize(const EquilibriumSpecs& specs) -> void
{
    pimpl->initialize(specs);
}

auto EquilibriumSensitivity::dndc(const String& cid) const -> VectorXdConstRef
{
    const auto idx = pimpl->params.index(cid);
    return pimpl->dndc.col(idx);
}

auto EquilibriumSensitivity::dndc(const Param& param) const -> VectorXdConstRef
{
    return dndc(param.id());
}

auto EquilibriumSensitivity::dndc() const -> MatrixXdConstRef
{
    return pimpl->dndc;
}

auto EquilibriumSensitivity::dndc(MatrixXdConstRef data) -> void
{
    auto& dndc = pimpl->dndc;
    errorif(dndc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dndc(MatrixXdConstRef).");
    errorif(dndc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dndc(MatrixXdConstRef).");
    dndc = data;
}

auto EquilibriumSensitivity::dpdc(const String& cid) const -> VectorXdConstRef
{
    const auto idx = pimpl->params.index(cid);
    return pimpl->dpdc.col(idx);
}

auto EquilibriumSensitivity::dpdc(const Param& param) const -> VectorXdConstRef
{
    return dpdc(param.id());
}

auto EquilibriumSensitivity::dpdc() const -> MatrixXdConstRef
{
    return pimpl->dpdc;
}

auto EquilibriumSensitivity::dpdc(MatrixXdConstRef data) -> void
{
    auto& dpdc = pimpl->dpdc;
    errorif(dpdc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dpdc(MatrixXdConstRef).");
    errorif(dpdc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dpdc(MatrixXdConstRef).");
    dpdc = data;
}

auto EquilibriumSensitivity::dqdc(const String& cid) const -> VectorXdConstRef
{
    const auto idx = pimpl->params.index(cid);
    return pimpl->dqdc.col(idx);
}

auto EquilibriumSensitivity::dqdc(const Param& param) const -> VectorXdConstRef
{
    return dqdc(param.id());
}

auto EquilibriumSensitivity::dqdc() const -> MatrixXdConstRef
{
    return pimpl->dqdc;
}

auto EquilibriumSensitivity::dqdc(MatrixXdConstRef data) -> void
{
    auto& dqdc = pimpl->dqdc;
    errorif(dqdc.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dqdc(MatrixXdConstRef).");
    errorif(dqdc.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dqdc(MatrixXdConstRef).");
    dqdc = data;
}

auto EquilibriumSensitivity::dndb() const -> MatrixXdConstRef
{
    return pimpl->dndb;
}

auto EquilibriumSensitivity::dndb(MatrixXdConstRef data) -> void
{
    auto& dndb = pimpl->dndb;
    errorif(dndb.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dndb(MatrixXdConstRef).");
    errorif(dndb.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dndb(MatrixXdConstRef).");
    dndb = data;
}

auto EquilibriumSensitivity::dpdb() const -> MatrixXdConstRef
{
    return pimpl->dpdb;
}

auto EquilibriumSensitivity::dpdb(MatrixXdConstRef data) -> void
{
    auto& dpdb = pimpl->dpdb;
    errorif(dpdb.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dpdb(MatrixXdConstRef).");
    errorif(dpdb.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dpdb(MatrixXdConstRef).");
    dpdb = data;
}

auto EquilibriumSensitivity::dqdb() const -> MatrixXdConstRef
{
    return pimpl->dqdb;
}

auto EquilibriumSensitivity::dqdb(MatrixXdConstRef data) -> void
{
    auto& dqdb = pimpl->dqdb;
    errorif(dqdb.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dqdb(MatrixXdConstRef).");
    errorif(dqdb.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dqdb(MatrixXdConstRef).");
    dqdb = data;
}

} // namespace Reaktoro
