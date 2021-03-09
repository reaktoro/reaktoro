// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include "EquilibriumSensitivity.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

struct EquilibriumSensitivity::Impl
{
    /// The chemical system associated with the sensitivity derivatives.
    ChemicalSystem system;

    /// The input variables *w* in the chemical equilibrium problem specifications.
    Strings inputs;

    /// The derivatives of the species amounts *n* with respect to input variables *w*.
    MatrixXd dndw;

    /// The derivatives of the control variables *p* with respect to input variables *w*.
    MatrixXd dpdw;

    /// The derivatives of the control variables *q* with respect to input variables *w*.
    MatrixXd dqdw;

    /// The derivatives of the species amounts *n* with respect to component amounts *b*.
    MatrixXd dndb;

    /// The derivatives of the control variables *p* with respect to component amounts *b*.
    MatrixXd dpdb;

    /// The derivatives of the control variables *q* with respect to component amounts *b*.
    MatrixXd dqdb;

    /// The total derivatives of the chemical properties *u* with respect to input variables *w*.
    MatrixXd dudw;

    /// The total derivatives of the chemical properties *u* with respect to component amounts *b*.
    MatrixXd dudb;

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
        inputs = specs.inputs();

        const EquilibriumDims dims(specs);

        const auto Nw = dims.Nw;
        const auto Nn = dims.Nn;
        const auto Np = dims.Np;
        const auto Nq = dims.Nq;
        const auto Nb = dims.Nb;

        dndw.resize(Nn, Nw);
        dpdw.resize(Np, Nw);
        dqdw.resize(Nq, Nw);
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

auto EquilibriumSensitivity::dndw(const String& wid) const -> VectorXdConstRef
{
    const auto idx = index(pimpl->inputs, wid);
    return pimpl->dndw.col(idx);
}

auto EquilibriumSensitivity::dndw(const Param& param) const -> VectorXdConstRef
{
    return dndw(param.id());
}

auto EquilibriumSensitivity::dndw() const -> MatrixXdConstRef
{
    return pimpl->dndw;
}

auto EquilibriumSensitivity::dndw(MatrixXdConstRef data) -> void
{
    auto& dndw = pimpl->dndw;
    errorif(dndw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dndw(MatrixXdConstRef).");
    errorif(dndw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dndw(MatrixXdConstRef).");
    dndw = data;
}

auto EquilibriumSensitivity::dpdw(const String& wid) const -> VectorXdConstRef
{
    const auto idx = index(pimpl->inputs, wid);
    return pimpl->dpdw.col(idx);
}

auto EquilibriumSensitivity::dpdw(const Param& param) const -> VectorXdConstRef
{
    return dpdw(param.id());
}

auto EquilibriumSensitivity::dpdw() const -> MatrixXdConstRef
{
    return pimpl->dpdw;
}

auto EquilibriumSensitivity::dpdw(MatrixXdConstRef data) -> void
{
    auto& dpdw = pimpl->dpdw;
    errorif(dpdw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dpdw(MatrixXdConstRef).");
    errorif(dpdw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dpdw(MatrixXdConstRef).");
    dpdw = data;
}

auto EquilibriumSensitivity::dqdw(const String& wid) const -> VectorXdConstRef
{
    const auto idx = index(pimpl->inputs, wid);
    return pimpl->dqdw.col(idx);
}

auto EquilibriumSensitivity::dqdw(const Param& param) const -> VectorXdConstRef
{
    return dqdw(param.id());
}

auto EquilibriumSensitivity::dqdw() const -> MatrixXdConstRef
{
    return pimpl->dqdw;
}

auto EquilibriumSensitivity::dqdw(MatrixXdConstRef data) -> void
{
    auto& dqdw = pimpl->dqdw;
    errorif(dqdw.rows() != data.rows(), "Mismatch number of rows in call to EquilibriumSensitivity::dqdw(MatrixXdConstRef).");
    errorif(dqdw.cols() != data.cols(), "Mismatch number of cols in call to EquilibriumSensitivity::dqdw(MatrixXdConstRef).");
    dqdw = data;
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

auto EquilibriumSensitivity::dudw() const -> MatrixXdConstRef
{
    return pimpl->dudw;
}

auto EquilibriumSensitivity::dudb() const -> MatrixXdConstRef
{
    return pimpl->dudb;
}

auto EquilibriumSensitivity::dudw(MatrixXdConstRef data) -> void
{
    pimpl->dudw = data;
}

auto EquilibriumSensitivity::dudb(MatrixXdConstRef data) -> void
{
    pimpl->dudb = data;
}

} // namespace Reaktoro
