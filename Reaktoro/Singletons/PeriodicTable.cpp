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

#include "PeriodicTable.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>

namespace Reaktoro {
namespace detail {

/// The default elements for the PeriodicTable object.
const std::vector<Element> default_elements =
{
    Element({ "H"  , "Hydrogen"      , 1   , 0.001007940 , 2.20 }),
    Element({ "He" , "Helium"        , 2   , 0.004002602 , 0.00 }),
    Element({ "Li" , "Lithium"       , 3   , 0.006941000 , 0.98 }),
    Element({ "Be" , "Beryllium"     , 4   , 0.009012180 , 1.57 }),
    Element({ "B"  , "Boron"         , 5   , 0.010811000 , 2.04 }),
    Element({ "C"  , "Carbon"        , 6   , 0.012011000 , 2.55 }),
    Element({ "N"  , "Nitrogen"      , 7   , 0.014006740 , 3.04 }),
    Element({ "O"  , "Oxygen"        , 8   , 0.015999400 , 3.44 }),
    Element({ "F"  , "Fluorine"      , 9   , 0.018998403 , 3.98 }),
    Element({ "Ne" , "Neon"          , 10  , 0.020179700 , 0.00 }),
    Element({ "Na" , "Sodium"        , 11  , 0.022989768 , 0.93 }),
    Element({ "Mg" , "Magnesium"     , 12  , 0.024305000 , 1.31 }),
    Element({ "Al" , "Aluminum"      , 13  , 0.026981539 , 1.61 }),
    Element({ "Si" , "Silicon"       , 14  , 0.028085500 , 1.90 }),
    Element({ "P"  , "Phosphorus"    , 15  , 0.030973762 , 2.19 }),
    Element({ "S"  , "Sulfur"        , 16  , 0.032066000 , 2.58 }),
    Element({ "Cl" , "Chlorine"      , 17  , 0.035452700 , 3.16 }),
    Element({ "Ar" , "Argon"         , 18  , 0.039948000 , 0.00 }),
    Element({ "K"  , "Potassium"     , 19  , 0.039098300 , 0.82 }),
    Element({ "Ca" , "Calcium"       , 20  , 0.040078000 , 1.00 }),
    Element({ "Sc" , "Scandium"      , 21  , 0.044955910 , 1.36 }),
    Element({ "Ti" , "Titanium"      , 22  , 0.047880000 , 1.54 }),
    Element({ "V"  , "Vanadium"      , 23  , 0.050941500 , 1.63 }),
    Element({ "Cr" , "Chromium"      , 24  , 0.051996100 , 1.66 }),
    Element({ "Mn" , "Manganese"     , 25  , 0.054938050 , 1.55 }),
    Element({ "Fe" , "Iron"          , 26  , 0.055847000 , 1.83 }),
    Element({ "Co" , "Cobalt"        , 27  , 0.058933200 , 1.88 }),
    Element({ "Ni" , "Nickel"        , 28  , 0.058693400 , 1.91 }),
    Element({ "Cu" , "Copper"        , 29  , 0.063546000 , 1.90 }),
    Element({ "Zn" , "Zinc"          , 30  , 0.065390000 , 1.65 }),
    Element({ "Ga" , "Gallium"       , 31  , 0.069723000 , 1.81 }),
    Element({ "Ge" , "Germanium"     , 32  , 0.072610000 , 2.01 }),
    Element({ "As" , "Arsenic"       , 33  , 0.074921590 , 2.18 }),
    Element({ "Se" , "Selenium"      , 34  , 0.078960000 , 2.55 }),
    Element({ "Br" , "Bromine"       , 35  , 0.079904000 , 2.96 }),
    Element({ "Kr" , "Krypton"       , 36  , 0.083800000 , 0.00 }),
    Element({ "Rb" , "Rubidium"      , 37  , 0.085467800 , 0.82 }),
    Element({ "Sr" , "Strontium"     , 38  , 0.087620000 , 0.95 }),
    Element({ "Y"  , "Yttrium"       , 39  , 0.088905850 , 1.22 }),
    Element({ "Zr" , "Zirconium"     , 40  , 0.091224000 , 1.33 }),
    Element({ "Nb" , "Niobium"       , 41  , 0.092906380 , 1.60 }),
    Element({ "Mo" , "Molybdenum"    , 42  , 0.095940000 , 2.16 }),
    Element({ "Tc" , "Technetium"    , 43  , 0.097907200 , 1.90 }),
    Element({ "Ru" , "Ruthenium"     , 44  , 0.101070000 , 2.20 }),
    Element({ "Rh" , "Rhodium"       , 45  , 0.102905500 , 2.28 }),
    Element({ "Pd" , "Palladium"     , 46  , 0.106420000 , 2.20 }),
    Element({ "Ag" , "Silver"        , 47  , 0.107868200 , 1.93 }),
    Element({ "Cd" , "Cadmium"       , 48  , 0.112411000 , 1.69 }),
    Element({ "In" , "Indium"        , 49  , 0.114818000 , 1.78 }),
    Element({ "Sn" , "Tin"           , 50  , 0.118710000 , 1.96 }),
    Element({ "Sb" , "Antimony"      , 51  , 0.121760000 , 2.05 }),
    Element({ "Te" , "Tellurium"     , 52  , 0.127600000 , 2.10 }),
    Element({ "I"  , "Iodine"        , 53  , 0.126904470 , 2.66 }),
    Element({ "Xe" , "Xenon"         , 54  , 0.131290000 , 0.00 }),
    Element({ "Cs" , "Cesium"        , 55  , 0.132905430 , 0.79 }),
    Element({ "Ba" , "Barium"        , 56  , 0.137327000 , 0.89 }),
    Element({ "La" , "Lanthanum"     , 57  , 0.138905500 , 1.10 }),
    Element({ "Ce" , "Cerium"        , 58  , 0.140115000 , 1.12 }),
    Element({ "Pr" , "Praseodymium"  , 59  , 0.140907650 , 1.13 }),
    Element({ "Nd" , "Neodymium"     , 60  , 0.144240000 , 1.14 }),
    Element({ "Pm" , "Promethium"    , 61  , 0.144912700 , 0.00 }),
    Element({ "Sm" , "Samarium"      , 62  , 0.150360000 , 1.17 }),
    Element({ "Eu" , "Europium"      , 63  , 0.151965000 , 0.00 }),
    Element({ "Gd" , "Gadolinium"    , 64  , 0.157250000 , 1.20 }),
    Element({ "Tb" , "Terbium"       , 65  , 0.158925340 , 1.20 }),
    Element({ "Dy" , "Dysprosium"    , 66  , 0.162500000 , 0.00 }),
    Element({ "Ho" , "Holmium"       , 67  , 0.164930320 , 1.23 }),
    Element({ "Er" , "Erbium"        , 68  , 0.167260000 , 1.24 }),
    Element({ "Tm" , "Thulium"       , 69  , 0.168934210 , 1.25 }),
    Element({ "Yb" , "Ytterbium"     , 70  , 0.173040000 , 1.10 }),
    Element({ "Lu" , "Lutetium"      , 71  , 0.174967000 , 1.27 }),
    Element({ "Hf" , "Hafnium"       , 72  , 0.178490000 , 1.30 }),
    Element({ "Ta" , "Tantalum"      , 73  , 0.180947900 , 1.50 }),
    Element({ "W"  , "Tungsten"      , 74  , 0.183840000 , 1.70 }),
    Element({ "Re" , "Rhenium"       , 75  , 0.186207000 , 1.90 }),
    Element({ "Os" , "Osmium"        , 76  , 0.190230000 , 2.20 }),
    Element({ "Ir" , "Iridium"       , 77  , 0.192220000 , 2.20 }),
    Element({ "Pt" , "Platinum"      , 78  , 0.195080000 , 2.28 }),
    Element({ "Au" , "Gold"          , 79  , 0.196966540 , 2.54 }),
    Element({ "Hg" , "Mercury"       , 80  , 0.200590000 , 2.00 }),
    Element({ "Tl" , "Thallium"      , 81  , 0.204383300 , 1.62 }),
    Element({ "Pb" , "Lead"          , 82  , 0.207200000 , 1.80 }),
    Element({ "Bi" , "Bismuth"       , 83  , 0.208980370 , 2.02 }),
    Element({ "Po" , "Polonium"      , 84  , 0.208982400 , 2.00 }),
    Element({ "At" , "Astatine"      , 85  , 0.209987100 , 2.20 }),
    Element({ "Rn" , "Radon"         , 86  , 0.222017600 , 0.00 }),
    Element({ "Fr" , "Francium"      , 87  , 0.223019700 , 0.70 }),
    Element({ "Ra" , "Radium"        , 88  , 0.226025400 , 0.90 }),
    Element({ "Ac" , "Actinium"      , 89  , 0.227027800 , 1.10 }),
    Element({ "Th" , "Thorium"       , 90  , 0.232038100 , 1.30 }),
    Element({ "Pa" , "Protactinium"  , 91  , 0.231035880 , 1.50 }),
    Element({ "U"  , "Uranium"       , 92  , 0.238028900 , 1.38 }),
    Element({ "Np" , "Neptunium"     , 93  , 0.237048000 , 1.36 }),
    Element({ "Pu" , "Plutonium"     , 94  , 0.244064200 , 1.28 }),
    Element({ "Am" , "Americium"     , 95  , 0.243061400 , 1.30 }),
    Element({ "Cm" , "Curium"        , 96  , 0.247070300 , 1.30 }),
    Element({ "Bk" , "Berkelium"     , 97  , 0.247070300 , 1.30 }),
    Element({ "Cf" , "Californium"   , 98  , 0.251079600 , 1.30 }),
    Element({ "Es" , "Einsteinium"   , 99  , 0.252083000 , 1.30 }),
    Element({ "Fm" , "Fermium"       , 100 , 0.257095100 , 1.30 }),
    Element({ "Md" , "Mendelevium"   , 101 , 0.258100000 , 1.30 }),
    Element({ "No" , "Nobelium"      , 102 , 0.259100900 , 1.30 }),
    Element({ "Lr" , "Lawrencium"    , 103 , 0.262110000 , 0.00 }),
    Element({ "Rf" , "Rutherfordium" , 104 , 0.261000000 , 0.00 }),
    Element({ "Db" , "Dubnium"       , 105 , 0.262000000 , 0.00 }),
    Element({ "Sg" , "Seaborgium"    , 106 , 0.266000000 , 0.00 }),
    Element({ "Bh" , "Bohrium"       , 107 , 0.264000000 , 0.00 }),
    Element({ "Hs" , "Hassium"       , 108 , 0.269000000 , 0.00 }),
    Element({ "Mt" , "Meitnerium"    , 109 , 0.268000000 , 0.00 }),
    Element({ "Ds" , "Darmstadtium"  , 110 , 0.269000000 , 0.00 }),
    Element({ "Rg" , "Roentgenium"   , 111 , 0.272000000 , 0.00 }),
    Element({ "Cn" , "Copernicium"   , 112 , 0.277000000 , 0.00 }),
    Element({ "Nh" , "Nihonium"      , 113 , 0.000000000 , 0.00 }),
    Element({ "Fl" , "Flerovium"     , 114 , 0.289000000 , 0.00 }),
    Element({ "Mc" , "Moscovium"     , 115 , 0.000000000 , 0.00 }),
    Element({ "Lv" , "Livermorium"   , 116 , 0.000000000 , 0.00 }),
    Element({ "Ts" , "Tennessine"    , 117 , 0.000000000 , 0.00 }),
    Element({ "Og" , "Oganesson"     , 118 , 0.000000000 , 0.00 }),
};

} // namespace detail

PeriodicTable::PeriodicTable()
: m_elements(detail::default_elements)
{}

PeriodicTable::~PeriodicTable()
{}

auto PeriodicTable::instance() -> PeriodicTable&
{
    static PeriodicTable obj;
    return obj;
}

auto PeriodicTable::elements() -> const std::vector<Element>&
{
    return instance().m_elements;
}

auto PeriodicTable::append(Element element) -> void
{
    auto& elements = instance().m_elements;
    elements.emplace_back(std::move(element));
}

auto PeriodicTable::size() -> std::size_t
{
    return elements().size();
}

auto PeriodicTable::elementWithName(std::string name) -> std::optional<Element>
{
    const auto idx = indexfn(elements(), [&](auto&& e) { return e.name() == name; });
    if(idx < size()) return elements()[idx];
    return {};
}

auto PeriodicTable::elementWithSymbol(std::string symbol) -> std::optional<Element>
{
    const auto idx = indexfn(elements(), [&](auto&& e) { return e.symbol() == symbol; });
    if(idx < size()) return elements()[idx];
    return {};
}

auto PeriodicTable::elementsWithTag(std::string tag) -> std::vector<Element>
{
    return filter(elements(), [&](auto&& e) { return contains(e.tags(), tag); });
}

auto PeriodicTable::elementsWithTags(const StringList& tags) -> std::vector<Element>
{
    return filter(elements(), [&](auto&& e) { return contained(tags, e.tags()); });
}

auto PeriodicTable::begin() const
{
    return m_elements.begin();
}

auto PeriodicTable::begin()
{
    return m_elements.begin();
}

auto PeriodicTable::end() const
{
    return m_elements.end();
}

auto PeriodicTable::end()
{
    return m_elements.end();
}

} // namespace Reaktoro
