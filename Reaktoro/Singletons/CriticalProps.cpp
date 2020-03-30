// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "CriticalProps.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {
namespace detail {

const std::deque<SubstanceCriticalProps> default_substances =
{
//   name                formula    Tcr       Pcr         omega
    {"Argon"           , "Ar"   , { 150.90 ,  48.98e+5 ,  0.0000 }},
    {"Ethylene"        , "C2H4" , { 282.30 ,  50.40e+5 ,  0.0870 }},
    {"Phenol"          , "C6H6O", { 694.30 ,  61.30e+5 ,  0.4440 }},
    {"Meta-Cresol"     , "C7H8O", {   0.00 ,   0.00e+5 ,  0.0000 }},
    {"Ortho-Cresol"    , "C7H8O", {   0.00 ,   0.00e+5 ,  0.0000 }},
    {"Para-Cresol"     , "C7H8O", {   0.00 ,   0.00e+5 ,  0.0000 }},
    {"Methane"         , "CH4"  , { 190.60 ,  45.99e+5 ,  0.0120 }},
    {"Carbon-Monoxide" , "CO"   , { 132.90 ,  34.99e+5 ,  0.0480 }},
    {"Carbon-Dioxide"  , "CO2"  , { 304.20 ,  73.83e+5 ,  0.2240 }},
    {"Hydrogen"        , "H2"   , {  33.19 ,  13.13e+5 , -0.2160 }},
    {"Steam"           , "H2O"  , { 647.10 , 220.55e+5 ,  0.3450 }},
    {"Hydrogen-Sulfide", "H2S"  , { 373.50 ,  89.63e+5 ,  0.0940 }},
    {"Helium"          , "He"   , {   5.20 ,   2.28e+5 , -0.3900 }},
    {"Krypton"         , "Kr"   , { 209.40 ,  55.02e+5 ,  0.0000 }},
    {"Nitrogen"        , "N2"   , { 126.20 ,  34.00e+5 ,  0.0380 }},
    {"Nitrous-Oxide"   , "N2O"  , { 309.60 ,  72.45e+5 ,  0.1410 }},
    {"Neon"            , "Ne"   , {  44.00 ,  27.00e+5 ,  0.0000 }},
    {"Ammonia"         , "NH3"  , { 405.70 , 112.80e+5 ,  0.2530 }},
    {"Nitric-Oxide"    , "NO"   , { 180.20 ,  64.80e+5 ,  0.5830 }},
    {"Oxygen"          , "O2"   , { 154.60 ,  50.43e+5 ,  0.0220 }},
    {"Radon"           , "Rn"   , { 377.00 ,  62.80e+5 ,  0.0000 }},
    {"Sulfur"          , "S2"   , {   0.00 ,   0.00e+5 ,  0.0000 }},
    {"Sulfur-Dioxide"  , "SO2"  , { 430.80 ,  78.84e+5 ,  0.2450 }},
    {"Xenon"           , "Xe"   , { 289.70 ,  58.40e+5 ,  0.0000 }},
};

} // namespace detail

SubstanceCriticalProps::SubstanceCriticalProps()
{
}

SubstanceCriticalProps::SubstanceCriticalProps(std::string name, const ChemicalFormula& formula)
: SubstanceCriticalProps(name, formula, {})
{
}

SubstanceCriticalProps::SubstanceCriticalProps(std::string name, const ChemicalFormula& formula, const SubstanceCriticalPropsData& data)
: m_name(name), m_formula(formula), m_data(data)
{
}

auto SubstanceCriticalProps::setTemperature(real value) -> void
{
    error(value <= 0.0, "Cannot set non-positive critical temperature value (", value, " K) to substance ", m_name, ".");
    m_data.Tcr = value;
}

auto SubstanceCriticalProps::setTemperature(real value, std::string unit) -> void
{
    setTemperature(units::convert(value, unit, "K"));
}

auto SubstanceCriticalProps::setPressure(real value) -> void
{
    error(value <= 0.0, "Cannot set non-positive critical pressure value (", value, " Pa) to substance ", m_name, ".");
    m_data.Pcr = value;
}

auto SubstanceCriticalProps::setPressure(real value, std::string unit) -> void
{
    error(value <= 0.0, "Cannot set non-positive critical pressure value (", value, " ", unit, ") to substance ", m_name, ".");
    m_data.Pcr = units::convert(value, unit, "Pa");
}

auto SubstanceCriticalProps::setAcentricFactor(real value) -> void
{
    m_data.omega = value;
}

auto SubstanceCriticalProps::name() const -> const std::string&
{
    return m_name;
}

auto SubstanceCriticalProps::formula() const -> const ChemicalFormula&
{
    return m_formula;
}

auto SubstanceCriticalProps::temperature() const -> real
{
    return m_data.Tcr;
}

auto SubstanceCriticalProps::pressure() const -> real
{
    return m_data.Pcr;
}

auto SubstanceCriticalProps::acentricFactor() const -> real
{
    return m_data.omega;
}

auto SubstanceCriticalProps::data() const -> const SubstanceCriticalPropsData&
{
    return m_data;
}

SubstanceCriticalProps::operator SubstanceCriticalPropsData() const
{
    return data();
}

CriticalProps::CriticalProps()
: m_substances(detail::default_substances)
{}

CriticalProps::~CriticalProps()
{}

auto CriticalProps::instance() -> CriticalProps&
{
    static CriticalProps obj;
    return obj;
}

auto CriticalProps::substances() -> const std::deque<SubstanceCriticalProps>&
{
    return instance().m_substances;
}

auto CriticalProps::append(SubstanceCriticalProps substance) -> void
{
    // Ensure there are no equivalent substances in the database (same name or equivalent chemical formulas).
    auto& substances = instance().m_substances;
    for(auto& current : substances)
    {
        if(substance.name == current.name || substance.formula.equivalent(current.formula)) {
            current = substance;
            return;
        }
    }
    substances.push_back(substance);
}

auto CriticalProps::size() -> std::size_t
{
    return substances().size();
}

auto CriticalProps::getWithName(const std::string& name) -> std::optional<SubstanceCriticalProps>
{
    const auto idx = indexfn(substances(), [&](auto&& s) { return s.name == name; });
    if(idx < size()) return substances()[idx];
    return {};
}

auto CriticalProps::getWithFormula(const ChemicalFormula& formula) -> std::optional<SubstanceCriticalProps>
{
    const auto idx = indexfn(substances(), [&](auto&& s) { return s.formula.equivalent(formula); });
    if(idx < size()) return substances()[idx];
    return {};
}

auto CriticalProps::get(const std::string& name_or_formula) -> std::optional<SubstanceCriticalProps>
{
    const auto subs = getWithName(name_or_formula);
    if(subs.has_value()) return subs;
    return getWithFormula(name_or_formula);
}

auto CriticalProps::begin() const
{
    return m_substances.begin();
}

auto CriticalProps::begin()
{
    return m_substances.begin();
}

auto CriticalProps::end() const
{
    return m_substances.end();
}

auto CriticalProps::end()
{
    return m_substances.end();
}

} // namespace Reaktoro