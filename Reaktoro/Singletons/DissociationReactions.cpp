// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "DissociationReactions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>

namespace Reaktoro {
namespace detail {

const Deque<DissociationReaction> default_reactions =
{
    { "AgCl", {{1, "Ag+"}, {1, "Cl-"}} },
    { "AgF", {{1, "Ag+"}, {1, "F-"}} },
    { "AgNO3", {{1, "Ag+"}, {1, "NO3-"}} },
    { "AlF3", {{1, "Al+++"}, {3, "F-"}} },
    { "As(OH)3", {{1, "As+++"}, {3, "OH-"}} },
    { "CaCl2", {{1, "Ca++"}, {2, "Cl-"}} },
    { "CaCO3", {{1, "Ca++"}, {1, "CO3--"}} },
    { "CaHPO4", {{1, "Ca++"}, {1, "HPO4--"}} },
    { "CaSO4", {{1, "Ca++"}, {1, "SO4--"}} },
    { "Cd(CN)2", {{1, "Cd++"}, {2, "CN-"}} },
    { "Cd(N3)2", {{1, "Cd++"}, {2, "N3-"}} },
    { "Cd(SCN)2", {{1, "Cd++"}, {2, "SCN-"}} },
    { "CdBr2", {{1, "Cd++"}, {2, "Br-"}} },
    { "CdCl2", {{1, "Cd++"}, {2, "Cl-"}} },
    { "CdI2", {{1, "Cd++"}, {2, "I-"}} },
    { "CdSeO4", {{1, "Cd++"}, {1, "SeO4--"}} },
    { "CdSO4", {{1, "Cd++"}, {1, "SO4--"}} },
    { "CeCl3", {{1, "Ce+++"}, {3, "Cl-"}} },
    { "CeF3", {{1, "Ce+++"}, {3, "F-"}} },
    { "Co(HS)2", {{1, "Co++"}, {2, "HS-"}} },
    { "CoBr2", {{1, "Co++"}, {2, "Br-"}} },
    { "CoI2", {{1, "Co++"}, {2, "I-"}} },
    { "CoSeO4", {{1, "Co++"}, {1, "SeO4--"}} },
    { "CoSO4", {{1, "Co++"}, {1, "SO4--"}} },
    { "CsBr", {{1, "Cs+"}, {1, "Br-"}} },
    { "CsCl", {{1, "Cs+"}, {1, "Cl-"}} },
    { "CsI", {{1, "Cs+"}, {1, "I-"}} },
    { "Cu(NO2)2", {{1, "Cu++"}, {2, "NO2-"}} },
    { "CuCl2", {{1, "Cu++"}, {2, "Cl-"}} },
    { "CuHPO4", {{1, "Cu++"}, {1, "HPO4--"}} },
    { "CuSO4", {{1, "Cu++"}, {1, "SO4--"}} },
    { "DyCl3", {{1, "Dy+++"}, {3, "Cl-"}} },
    { "DyF3", {{1, "Dy+++"}, {3, "F-"}} },
    { "ErCl3", {{1, "Er+++"}, {3, "Cl-"}} },
    { "ErF3", {{1, "Er+++"}, {3, "F-"}} },
    { "EuCl2", {{1, "Eu++"}, {2, "Cl-"}} },
    { "EuCl3", {{1, "Eu+++"}, {3, "Cl-"}} },
    { "EuF2", {{1, "Eu++"}, {2, "F-"}} },
    { "EuF3", {{1, "Eu+++"}, {3, "F-"}} },
    { "FeCl2", {{1, "Fe++"}, {2, "Cl-"}} },
    { "FeHPO4", {{1, "Fe++"}, {1, "HPO4--"}} },
    { "FeSO4", {{1, "Fe++"}, {1, "SO4--"}} },
    { "GdCl3", {{1, "Gd+++"}, {3, "Cl-"}} },
    { "GdF3", {{1, "Gd+++"}, {3, "F-"}} },
    { "H2CrO4", {{2, "H+"}, {1, "CrO4--"}} },
    { "H2F2", {{2, "H+"}, {2, "F-"}} },
    { "H2S", {{1, "H+"}, {1, "HS-"}} },
    { "H2Se", {{2, "H+"}, {1, "Se--"}} },
    { "H2SeO3", {{2, "H+"}, {1, "SeO3--"}} },
    { "H2SO3", {{2, "H+"}, {1, "SO3--"}} },
    { "H2SO4", {{2, "H+"}, {1, "SO4--"}} },
    { "H2TcO4", {{2, "H+"}, {1, "TcO4--"}} },
    { "H3AsO4", {{1, "H+"}, {1, "H2AsO4-"}} },
    { "H3PO4", {{2, "H+"}, {1, "HPO4--"}} },
    { "HBrO", {{1, "H+"}, {1, "BrO-"}} },
    { "HCl", {{1, "H+"}, {1, "Cl-"}} },
    { "HClO", {{1, "H+"}, {1, "ClO-"}} },
    { "HClO2", {{1, "H+"}, {1, "ClO2-"}} },
    { "HCN", {{1, "H+"}, {1, "CN-"}} },
    { "HF", {{1, "H+"}, {1, "F-"}} },
    { "HIO3", {{1, "H+"}, {1, "IO3-"}} },
    { "HN3", {{1, "H+"}, {1, "N3-"}} },
    { "HNO2", {{1, "H+"}, {1, "NO2-"}} },
    { "HNO3", {{1, "H+"}, {1, "NO3-"}} },
    { "HoCl3", {{1, "Ho+++"}, {3, "Cl-"}} },
    { "HoF3", {{1, "Ho+++"}, {3, "F-"}} },
    { "KBr", {{1, "K+"}, {1, "Br-"}} },
    { "KCl", {{1, "K+"}, {1, "Cl-"}} },
    { "KHSO4", {{1, "K+"}, {1, "H+"}, {1, "SO4--"}} },
    { "KI", {{1, "K+"}, {1, "I-"}} },
    { "LaCl3", {{1, "La+++"}, {3, "Cl-"}} },
    { "LaF3", {{1, "La+++"}, {3, "F-"}} },
    { "LiCl", {{1, "Li+"}, {1, "Cl-"}} },
    { "LuCl3", {{1, "Lu+++"}, {3, "Cl-"}} },
    { "LuF3", {{1, "Lu+++"}, {3, "F-"}} },
    { "MgCl2", {{1, "Mg++"}, {2, "Cl-"}} },
    { "MgCO3", {{1, "Mg++"}, {1, "CO3--"}} },
    { "MgHPO4", {{1, "Mg++"}, {1, "HPO4--"}} },
    { "MgSO4", {{1, "Mg++"}, {1, "SO4--"}} },
    { "Mn(NO3)2", {{1, "Mn++"}, {2, "NO3-"}} },
    { "MnHPO4", {{1, "Mn++"}, {1, "HPO4--"}} },
    { "MnSeO4", {{1, "Mn++"}, {1, "SeO4--"}} },
    { "MnSO4", {{1, "Mn++"}, {1, "SO4--"}} },
    { "NaBr", {{1, "Na+"}, {1, "Br-"}} },
    { "NaCl", {{1, "Na+"}, {1, "Cl-"}} },
    { "NaF", {{1, "Na+"}, {1, "F-"}} },
    { "NaHCO3", {{1, "Na+"}, {1, "HCO3-"}} },
    { "NaI", {{1, "Na+"}, {1, "I-"}} },
    { "NaOH", {{1, "Na+"}, {1, "OH-"}} },
    { "NdCl3", {{1, "Nd+++"}, {3, "Cl-"}} },
    { "NdF3", {{1, "Nd+++"}, {3, "F-"}} },
    { "Ni(NO3)2", {{1, "Ni++"}, {2, "NO3-"}} },
    { "NiSeO4", {{1, "Ni++"}, {1, "SeO4--"}} },
    { "NiSO4", {{1, "Ni++"}, {1, "SO4--"}} },
    { "Np(H2PO4)3", {{1, "Np+++"}, {3, "H+"}, {3, "HPO4--"}} },
    { "Np(HPO4)2", {{1, "Np++++"}, {2, "HPO4--"}} },
    { "Np(SO4)2", {{1, "Np++++"}, {2, "SO4--"}} },
    { "NpO2Cl", {{1, "NpO2+"}, {1, "Cl-"}} },
    { "NpO2F", {{1, "NpO2+"}, {1, "F-"}} },
    { "NpO2F2", {{1, "NpO2++"}, {2, "F-"}} },
    { "NpO2H2PO4", {{1, "NpO2+"}, {1, "H+"}, {1, "HPO4--"}} },
    { "NpO2HPO4", {{1, "NpO2++"}, {1, "HPO4--"}} },
    { "NpO2SO4", {{1, "NpO2++"}, {1, "SO4--"}} },
    { "Pb(BrO3)2", {{1, "Pb++"}, {2, "BrO3-"}} },
    { "Pb(ClO3)2", {{1, "Pb++"}, {2, "ClO3-"}} },
    { "Pb(SCN)2", {{1, "Pb++"}, {2, "SCN-"}} },
    { "PbBr2", {{1, "Pb++"}, {2, "Br-"}} },
    { "PbCl2", {{1, "Pb++"}, {2, "Cl-"}} },
    { "PbF2", {{1, "Pb++"}, {2, "F-"}} },
    { "PbHPO4", {{1, "Pb++"}, {1, "HPO4--"}} },
    { "PbI2", {{1, "Pb++"}, {2, "I-"}} },
    { "PdCl2", {{1, "Pd++"}, {2, "Cl-"}} },
    { "PrCl3", {{1, "Pr+++"}, {3, "Cl-"}} },
    { "PrF3", {{1, "Pr+++"}, {3, "F-"}} },
    { "Pu(HPO4)2", {{1, "Pu++++"}, {2, "HPO4--"}} },
    { "Pu(SO4)2", {{1, "Pu++++"}, {2, "SO4--"}} },
    { "PuF4", {{1, "Pu++++"}, {4, "F-"}} },
    { "PuO2F2", {{1, "PuO2++"}, {2, "F-"}} },
    { "PuO2SO4", {{1, "PuO2++"}, {1, "SO4--"}} },
    { "RbBr", {{1, "Rb+"}, {1, "Br-"}} },
    { "RbCl", {{1, "Rb+"}, {1, "Cl-"}} },
    { "RbF", {{1, "Rb+"}, {1, "F-"}} },
    { "RbI", {{1, "Rb+"}, {1, "I-"}} },
    { "RuCl3", {{1, "Ru+++"}, {3, "Cl-"}} },
    { "Ru(OH)2Cl2", {{1, "Ru(OH)2++"}, {2, "Cl-"}} },
    { "Ru(OH)2SO4", {{1, "Ru(OH)2++"}, {1, "SO4--"}} },
    { "RuSO4", {{1, "Ru++"}, {1, "SO4--"}} },
    { "SmCl3", {{1, "Sm+++"}, {3, "Cl-"}} },
    { "SmF3", {{1, "Sm+++"}, {3, "F-"}} },
    { "Sn(SO4)2", {{1, "Sn++++"}, {2, "SO4--"}} },
    { "SnCl2", {{1, "Sn++"}, {2, "Cl-"}} },
    { "SnF2", {{1, "Sn++"}, {2, "F-"}} },
    { "SrHPO4", {{1, "Sr++"}, {1, "HPO4--"}} },
    { "SrSO4", {{1, "Sr++"}, {1, "SO4--"}} },
    { "TbCl3", {{1, "Tb+++"}, {3, "Cl-"}} },
    { "TbF3", {{1, "Tb+++"}, {3, "F-"}} },
    { "Th(HPO4)2", {{1, "Th++++"}, {2, "HPO4--"}} },
    { "Th(SO4)2", {{1, "Th++++"}, {2, "SO4--"}} },
    { "ThCl4", {{1, "Th++++"}, {4, "Cl-"}} },
    { "ThF4", {{1, "Th++++"}, {4, "F-"}} },
    { "TmCl3", {{1, "Tm+++"}, {3, "Cl-"}} },
    { "TmF3", {{1, "Tm+++"}, {3, "F-"}} },
    { "U(SO4)2", {{1, "U++++"}, {2, "SO4--"}} },
    { "UF4", {{1, "U++++"}, {4, "F-"}} },
    { "UO2(H2PO4)2", {{1, "UO2++"}, {2, "H+"}, {2, "HPO4--"}} },
    { "UO2(IO3)2", {{1, "UO2++"}, {2, "IO3-"}} },
    { "UO2(N3)2", {{1, "UO2++"}, {2, "N3-"}} },
    { "UO2(SCN)2", {{1, "UO2++"}, {2, "SCN-"}} },
    { "UO2Cl2", {{1, "UO2++"}, {2, "Cl-"}} },
    { "UO2F2", {{1, "UO2++"}, {2, "F-"}} },
    { "UO2HPO4", {{1, "UO2++"}, {1, "HPO4--"}} },
    { "UO2SO3", {{1, "UO2++"}, {1, "SO3--"}} },
    { "UO2SO4", {{1, "UO2++"}, {1, "SO4--"}} },
    { "VO2F", {{1, "VO2+"}, {1, "F-"}} },
    { "VO2H2PO4", {{1, "VO2+"}, {1, "H+"}, {1, "HPO4--"}} },
    { "VOF2", {{1, "VO++"}, {2, "F-"}} },
    { "VOSO4", {{1, "VO++"}, {1, "SO4--"}} },
    { "YbCl3", {{1, "Yb+++"}, {3, "Cl-"}} },
    { "YbF3", {{1, "Yb+++"}, {3, "F-"}} },
    { "YF3", {{1, "Y+++"}, {3, "F-"}} },
    { "Zn(N3)2", {{1, "Zn++"}, {2, "N3-"}} },
    { "Zn(SCN)2", {{1, "Zn++"}, {2, "SCN-"}} },
    { "ZnBr2", {{1, "Zn++"}, {2, "Br-"}} },
    { "ZnCl2", {{1, "Zn++"}, {2, "Cl-"}} },
    { "ZnHPO4", {{1, "Zn++"}, {1, "HPO4--"}} },
    { "ZnI2", {{1, "Zn++"}, {2, "I-"}} },
    { "ZnSeO4", {{1, "Zn++"}, {1, "SeO4--"}} },
    { "ZnSO4", {{1, "Zn++"}, {1, "SO4--"}} },
    { "Zr(SO4)2", {{1, "Zr++++"}, {2, "SO4--"}} },
    { "ZrF4", {{1, "Zr++++"}, {4, "F-"}} },
};

} // namespace detail

DissociationReactions::DissociationReactions()
: m_reactions(detail::default_reactions)
{}

DissociationReactions::~DissociationReactions()
{}

auto DissociationReactions::instance() -> DissociationReactions&
{
    static DissociationReactions obj;
    return obj;
}

auto DissociationReactions::reactions() -> const Deque<DissociationReaction>&
{
    return instance().m_reactions;
}

auto DissociationReactions::reset() -> void
{
    instance().m_reactions = detail::default_reactions;
}

auto DissociationReactions::append(DissociationReaction reaction) -> void
{
    // Ensure there are no equivalent complex substances in the database.
    auto& reactions = instance().m_reactions;
    for(auto& current : reactions)
    {
        if(reaction.complex.equivalent(current.complex)) {
            current.ions = reaction.ions;
            return;
        }
    }
    reactions.push_back(reaction);
}

auto DissociationReactions::size() -> std::size_t
{
    return reactions().size();
}

auto DissociationReactions::get(const ChemicalFormula& complex) -> std::optional<DissociationReaction>
{
    const auto idx = indexfn(reactions(), [&](auto&& r) { return r.complex.equivalent(complex); });
    if(idx < size()) return reactions()[idx];
    return {};
}

auto DissociationReactions::coefficient(const ChemicalFormula& complex, const ChemicalFormula& ion) -> double
{
    const auto reaction = get(complex);
    if(reaction.has_value())
        for(auto&& [coeff, product] : reaction.value().ions)
            if(product.equivalent(ion))
                return coeff;
    return 0.0;
}

auto DissociationReactions::begin() const
{
    return m_reactions.begin();
}

auto DissociationReactions::begin()
{
    return m_reactions.begin();
}

auto DissociationReactions::end() const
{
    return m_reactions.end();
}

auto DissociationReactions::end()
{
    return m_reactions.end();
}

} // namespace Reaktoro
