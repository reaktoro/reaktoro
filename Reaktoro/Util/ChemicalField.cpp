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
//
//#include "ChemicalField.hpp"
//
//// C++ includes
//#include <iomanip>
//
//// Reaktoro includes
//#include <Reaktoro/Core/ChemicalSystem.hpp>
//#include <Reaktoro/Core/Partition.hpp>
//#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
//
//namespace Reaktoro {
//
//struct ChemicalField::Impl
//{
//    /// The partition of the chemical system.
//    Partition partition;
//
//    /// The number of points in the field.
//    Index npoints = 0;
//
//    /// The values of the scalar chemical field.
//    VectorXr val;
//
//    /// The derivatives of the scalar chemical field with respect to temperature.
//    VectorXr ddT;
//
//    /// The derivatives of the scalar chemical field with respect to pressure.
//    VectorXr ddP;
//
//    /// The derivatives of the scalar chemical field with respect to the amounts of each equilibrium element.
//    std::vector<VectorXr> ddbe;
//
//    /// The derivatives of the scalar chemical field with respect to the amounts of each kinetic species.
//    std::vector<VectorXr> ddnk;
//
//    /// Auxiliary vectors to avoid recurrent memory allocation.
//    VectorXr scalar_ne, scalar_nk, scalar_be;
//
//    /// Construct a default ChemicalField instance.
//    Impl()
//    {}
//
//    /// Construct a ChemicalField instance with given chemical system partition.
//    Impl(const Partition& partition, Index npoints)
//    : partition(partition), npoints(npoints),
//      val(zeros(npoints)), ddT(zeros(npoints)), ddP(zeros(npoints)),
//      ddbe(partition.numEquilibriumElements(), zeros(npoints)),
//      ddnk(partition.numKineticSpecies(), zeros(npoints))
//    {}
//
//    /// Set the field at the i-th point with a real instance.
//    auto set(Index i, const real& scalar, const EquilibriumSensitivity& sensitivity) -> void
//    {
//        // The indices of the equilibrium and kinetic species
//        const Indices& ispecies_e = partition.indicesEquilibriumSpecies();
//        const Indices& ispecies_k = partition.indicesKineticSpecies();
//
//        // Auxiliary references to sensitivity values
//        VectorXrConstRef ne_T  = sensitivity.dnedT;
//        VectorXrConstRef ne_P  = sensitivity.dnedP;
//        MatrixXdConstRef ne_be = sensitivity.dnedbe;
//
//        // Extract the derivatives of scalar w.r.t. amounts of equilibrium species
//        scalar_ne = rows(scalar.ddn, ispecies_e);
//
//        // Extract the derivatives of scalar w.r.t. amounts of kinetic species
//        scalar_nk = rows(scalar.ddn, ispecies_k);
//
//        // Calculte the derivatives of scalar w.r.t. amounts of equilibrium elements
//        scalar_be = tr(scalar_ne) * ne_be;
//
//        // Set the i-th position of the scalar field with given scalar value
//        val[i] = scalar.val;
//
//        // Set derivative w.r.t. temperature at the i-th position
//        ddT[i] = scalar.ddT + dot(scalar_ne, ne_T);
//
//        // Set derivative w.r.t. pressure at the i-th position
//        ddP[i] = scalar.ddP + dot(scalar_ne, ne_P);
//
//        // Set derivative w.r.t. amounts of equilibrium elements at the i-th position
//        for(Index j = 0; j < ddbe.size(); ++j)
//            ddbe[j][i] = scalar_be[j];
//
//        // Set derivative w.r.t. amounts of kinetic species at the i-th position
//        for(Index j = 0; j < ddnk.size(); ++j)
//            ddnk[j][i] = scalar_nk[j];
//    }
//};
//
//ChemicalField::ChemicalField()
//: pimpl(new Impl())
//{}
//
//ChemicalField::ChemicalField(const Partition& partition, Index npoints)
//: pimpl(new Impl(partition, npoints))
//{}
//
//ChemicalField::ChemicalField(const ChemicalField& other)
//: pimpl(new Impl(*other.pimpl))
//{}
//
//ChemicalField::~ChemicalField()
//{}
//
//auto ChemicalField::operator=(ChemicalField other) -> ChemicalField&
//{
//    pimpl = std::move(other.pimpl);
//    return *this;
//}
//
//auto ChemicalField::set(Index i, const real& scalar, const EquilibriumSensitivity& sensitivity) -> void
//{
//    pimpl->set(i, scalar, sensitivity);
//}
//
//auto ChemicalField::partition() const -> const Partition&
//{
//    return pimpl->partition;
//}
//
//auto ChemicalField::size() const -> Index
//{
//    return pimpl->npoints;
//}
//
//auto ChemicalField::val() -> VectorXrRef
//{
//    return pimpl->val;
//}
//
//auto ChemicalField::val() const -> VectorXrConstRef
//{
//    return pimpl->val;
//}
//
//auto ChemicalField::ddT() -> VectorXrRef
//{
//    return pimpl->ddT;
//}
//
//auto ChemicalField::ddT() const -> VectorXrConstRef
//{
//    return pimpl->ddT;
//}
//
//auto ChemicalField::ddP() -> VectorXrRef
//{
//    return pimpl->ddP;
//}
//
//auto ChemicalField::ddP() const -> VectorXrConstRef
//{
//    return pimpl->ddP;
//}
//
//auto ChemicalField::ddbe() -> std::vector<VectorXr>&
//{
//    return pimpl->ddbe;
//}
//
//auto ChemicalField::ddbe() const -> const std::vector<VectorXr>&
//{
//    return pimpl->ddbe;
//}
//
//auto ChemicalField::ddnk() -> std::vector<VectorXr>&
//{
//    return pimpl->ddnk;
//}
//
//auto ChemicalField::ddnk() const -> const std::vector<VectorXr>&
//{
//    return pimpl->ddnk;
//}
//
//auto operator<<(std::ostream& out, const ChemicalField& f) -> std::ostream&
//{
//    const Partition& partition = f.partition();
//    const ChemicalSystem& system = partition.system();
//
//    const Indices& iee = partition.indicesEquilibriumElements();
//    const Indices& iks = partition.indicesKineticSpecies();
//
//    const Index Ee = partition.numEquilibriumElements();
//    const Index Nk = partition.numKineticSpecies();
//
//    out << std::left << std::setw(10) << "k";
//    out << std::left << std::setw(20) << "val";
//    out << std::left << std::setw(20) << "ddT";
//    out << std::left << std::setw(20) << "ddP";
//    for(Index i = 0; i < Ee; ++i)
//        out << std::left << std::setw(20) << "ddbe(" + system.element(iee[i]).name() + ")";
//    for(Index i = 0; i < Nk; ++i)
//        out << std::left << std::setw(20) << "ddnk(" + system.species(iks[i]).name() + ")";
//    out << std::endl;
//    for(Index k = 0; k < f.size(); ++k)
//    {
//        out << std::left << std::setw(10) << k;
//        out << std::left << std::setw(20) << f.val()[k];
//        out << std::left << std::setw(20) << f.ddT()[k];
//        out << std::left << std::setw(20) << f.ddP()[k];
//        for(Index i = 0; i < Ee; ++i)
//            out << std::left << std::setw(20) << f.ddbe()[i][k];
//        for(Index i = 0; i < Nk; ++i)
//            out << std::left << std::setw(20) << f.ddnk()[i][k];
//        out << std::endl;
//    }
//
//    return out;
//}
//
//}  // namespace Reaktoro
