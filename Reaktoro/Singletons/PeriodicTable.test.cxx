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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Singletons/PeriodicTable.hpp>
using namespace Reaktoro;

TEST_CASE("Testing PeriodicTable", "[PeriodicTable]")
{
    REQUIRE(PeriodicTable::size() == PeriodicTable::elements().size());

    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("H"  ).value().name() == "Hydrogen"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("He" ).value().name() == "Helium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Li" ).value().name() == "Lithium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Be" ).value().name() == "Beryllium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("B"  ).value().name() == "Boron"         );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("C"  ).value().name() == "Carbon"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("N"  ).value().name() == "Nitrogen"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("O"  ).value().name() == "Oxygen"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("F"  ).value().name() == "Fluorine"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ne" ).value().name() == "Neon"          );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Na" ).value().name() == "Sodium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Mg" ).value().name() == "Magnesium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Al" ).value().name() == "Aluminum"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Si" ).value().name() == "Silicon"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("P"  ).value().name() == "Phosphorus"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("S"  ).value().name() == "Sulfur"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cl" ).value().name() == "Chlorine"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ar" ).value().name() == "Argon"         );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("K"  ).value().name() == "Potassium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ca" ).value().name() == "Calcium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sc" ).value().name() == "Scandium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ti" ).value().name() == "Titanium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("V"  ).value().name() == "Vanadium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cr" ).value().name() == "Chromium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Mn" ).value().name() == "Manganese"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Fe" ).value().name() == "Iron"          );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Co" ).value().name() == "Cobalt"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ni" ).value().name() == "Nickel"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cu" ).value().name() == "Copper"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Zn" ).value().name() == "Zinc"          );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ga" ).value().name() == "Gallium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ge" ).value().name() == "Germanium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("As" ).value().name() == "Arsenic"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Se" ).value().name() == "Selenium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Br" ).value().name() == "Bromine"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Kr" ).value().name() == "Krypton"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Rb" ).value().name() == "Rubidium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sr" ).value().name() == "Strontium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Y"  ).value().name() == "Yttrium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Zr" ).value().name() == "Zirconium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Nb" ).value().name() == "Niobium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Mo" ).value().name() == "Molybdenum"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Tc" ).value().name() == "Technetium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ru" ).value().name() == "Ruthenium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Rh" ).value().name() == "Rhodium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pd" ).value().name() == "Palladium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ag" ).value().name() == "Silver"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cd" ).value().name() == "Cadmium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("In" ).value().name() == "Indium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sn" ).value().name() == "Tin"           );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sb" ).value().name() == "Antimony"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Te" ).value().name() == "Tellurium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("I"  ).value().name() == "Iodine"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Xe" ).value().name() == "Xenon"         );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cs" ).value().name() == "Cesium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ba" ).value().name() == "Barium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("La" ).value().name() == "Lanthanum"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ce" ).value().name() == "Cerium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pr" ).value().name() == "Praseodymium"  );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Nd" ).value().name() == "Neodymium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pm" ).value().name() == "Promethium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sm" ).value().name() == "Samarium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Eu" ).value().name() == "Europium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Gd" ).value().name() == "Gadolinium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Tb" ).value().name() == "Terbium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Dy" ).value().name() == "Dysprosium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ho" ).value().name() == "Holmium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Er" ).value().name() == "Erbium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Tm" ).value().name() == "Thulium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Yb" ).value().name() == "Ytterbium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Lu" ).value().name() == "Lutetium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Hf" ).value().name() == "Hafnium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ta" ).value().name() == "Tantalum"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("W"  ).value().name() == "Tungsten"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Re" ).value().name() == "Rhenium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Os" ).value().name() == "Osmium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ir" ).value().name() == "Iridium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pt" ).value().name() == "Platinum"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Au" ).value().name() == "Gold"          );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Hg" ).value().name() == "Mercury"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Tl" ).value().name() == "Thallium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pb" ).value().name() == "Lead"          );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Bi" ).value().name() == "Bismuth"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Po" ).value().name() == "Polonium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("At" ).value().name() == "Astatine"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Rn" ).value().name() == "Radon"         );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Fr" ).value().name() == "Francium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ra" ).value().name() == "Radium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ac" ).value().name() == "Actinium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Th" ).value().name() == "Thorium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pa" ).value().name() == "Protactinium"  );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("U"  ).value().name() == "Uranium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Np" ).value().name() == "Neptunium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Pu" ).value().name() == "Plutonium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Am" ).value().name() == "Americium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cm" ).value().name() == "Curium"        );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Bk" ).value().name() == "Berkelium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cf" ).value().name() == "Californium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Es" ).value().name() == "Einsteinium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Fm" ).value().name() == "Fermium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Md" ).value().name() == "Mendelevium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("No" ).value().name() == "Nobelium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Lr" ).value().name() == "Lawrencium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Rf" ).value().name() == "Rutherfordium" );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Db" ).value().name() == "Dubnium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Sg" ).value().name() == "Seaborgium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Bh" ).value().name() == "Bohrium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Hs" ).value().name() == "Hassium"       );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Mt" ).value().name() == "Meitnerium"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ds" ).value().name() == "Darmstadtium"  );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Rg" ).value().name() == "Roentgenium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Cn" ).value().name() == "Copernicium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Nh" ).value().name() == "Nihonium"      );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Fl" ).value().name() == "Flerovium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Mc" ).value().name() == "Moscovium"     );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Lv" ).value().name() == "Livermorium"   );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Ts" ).value().name() == "Tennessine"    );
    REQUIRE_NOTHROW( PeriodicTable::elementWithSymbol("Og" ).value().name() == "Oganesson"     );

    REQUIRE_NOTHROW(PeriodicTable::elementWithName("Hydrogen").value().symbol() == "H" );
    REQUIRE_NOTHROW(PeriodicTable::elementWithName("Helium"  ).value().symbol() == "He");
    REQUIRE_NOTHROW(PeriodicTable::elementWithName("Lithium" ).value().symbol() == "Li");

    REQUIRE_FALSE(PeriodicTable::elementWithSymbol("Ab"));
    REQUIRE_FALSE(PeriodicTable::elementWithSymbol("Xy"));

    PeriodicTable::append(Element().withSymbol("Aa").withTags({"tag1", "tag2", "tag3"}));
    PeriodicTable::append(Element().withSymbol("Bb").withTags({"tag1", "tag3"}));
    PeriodicTable::append(Element().withSymbol("Cc").withTags({"tag1"}));

    const auto elements_with_tag = PeriodicTable::elementsWithTag("tag1");

    REQUIRE(elements_with_tag.size() == 3);
    REQUIRE(elements_with_tag[0].symbol() == "Aa");
    REQUIRE(elements_with_tag[1].symbol() == "Bb");
    REQUIRE(elements_with_tag[2].symbol() == "Cc");

    const auto elements_with_tags = PeriodicTable::elementsWithTags({"tag3", "tag1"});

    REQUIRE(elements_with_tags.size() == 2);
    REQUIRE(elements_with_tags[0].symbol() == "Aa");
    REQUIRE(elements_with_tags[1].symbol() == "Bb");
}
