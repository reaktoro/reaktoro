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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Singletons/CriticalProps.hpp>
using namespace Reaktoro;

namespace withName {

auto get(std::string name) { return CriticalProps::get(name).value(); }
auto temperature(std::string name) { return get(name).temperature(); }
auto pressure(std::string name) { return get(name).pressure(); }
auto acentricFactor(std::string name) { return get(name).acentricFactor(); }

} // namespace withName

namespace withNames {

auto get(std::vector<std::string> names) { return CriticalProps::get(names).value(); }
auto temperature(std::vector<std::string> names) { return get(names).temperature(); }
auto pressure(std::vector<std::string> names) { return get(names).pressure(); }
auto acentricFactor(std::vector<std::string> names) { return get(names).acentricFactor(); }

} // namespace withNames

TEST_CASE("Testing CriticalProps", "[CriticalProps]")
{
    REQUIRE(CriticalProps::size() == CriticalProps::data().size());
    REQUIRE_NOTHROW( withName::get("methane")              );
    REQUIRE_NOTHROW( withName::get("ethane")               );
    REQUIRE_NOTHROW( withName::get("propane")              );
    REQUIRE_NOTHROW( withName::get("N-BUTANE")             );
    REQUIRE_NOTHROW( withName::get("N-PENTANE")            );
    REQUIRE_NOTHROW( withName::get("N-HEXANE")             );
    REQUIRE_NOTHROW( withName::get("N-HEPTANE")            );
    REQUIRE_NOTHROW( withName::get("N-OCTANE")             );
    REQUIRE_NOTHROW( withName::get("N-NONANE")             );
    REQUIRE_NOTHROW( withName::get("N-DECANE")             );
    REQUIRE_NOTHROW( withName::get("ISOBUTANE")            );
    REQUIRE_NOTHROW( withName::get("ISOOCTANE")            );
    REQUIRE_NOTHROW( withName::get("CYCLOPENTANE")         );
    REQUIRE_NOTHROW( withName::get("CYCLOHEXANE")          );
    REQUIRE_NOTHROW( withName::get("METHYLCYCLOPENTANE")   );
    REQUIRE_NOTHROW( withName::get("METHYLCYCLOHEXANE")    );
    REQUIRE_NOTHROW( withName::get("ETHYLENE")             );
    REQUIRE_NOTHROW( withName::get("PROPYLENE")            );
    REQUIRE_NOTHROW( withName::get("1-BUTENE")             );
    REQUIRE_NOTHROW( withName::get("CIS-2-BUTENE")         );
    REQUIRE_NOTHROW( withName::get("TRANS-2-BUTENE")       );
    REQUIRE_NOTHROW( withName::get("1-HEXENE")             );
    REQUIRE_NOTHROW( withName::get("ISOBUTYLENE")          );
    REQUIRE_NOTHROW( withName::get("1,3-BUTADIENE")        );
    REQUIRE_NOTHROW( withName::get("CYCLOHEXENE")          );
    REQUIRE_NOTHROW( withName::get("ACETYLENE")            );
    REQUIRE_NOTHROW( withName::get("BENZENE")              );
    REQUIRE_NOTHROW( withName::get("TOLUENE")              );
    REQUIRE_NOTHROW( withName::get("ETHYLBENZENE")         );
    REQUIRE_NOTHROW( withName::get("CUMENE")               );
    REQUIRE_NOTHROW( withName::get("O-XYLENE")             );
    REQUIRE_NOTHROW( withName::get("M-XYLENE")             );
    REQUIRE_NOTHROW( withName::get("P-XYLENE")             );
    REQUIRE_NOTHROW( withName::get("STYRENE")              );
    REQUIRE_NOTHROW( withName::get("NAPHTHALENE")          );
    REQUIRE_NOTHROW( withName::get("BIPHENYL")             );
    REQUIRE_NOTHROW( withName::get("FORMALDEHYDE")         );
    REQUIRE_NOTHROW( withName::get("ACETALDEHYDE")         );
    REQUIRE_NOTHROW( withName::get("METHYL-ACETATE")       );
    REQUIRE_NOTHROW( withName::get("ETHYL-ACETATE")        );
    REQUIRE_NOTHROW( withName::get("ACETONE")              );
    REQUIRE_NOTHROW( withName::get("METHYL-ETHYL-KETONE")  );
    REQUIRE_NOTHROW( withName::get("DIETHYL-ETHER")        );
    REQUIRE_NOTHROW( withName::get("METHYL-T-BUTYL-ETHER") );
    REQUIRE_NOTHROW( withName::get("METHANOL")             );
    REQUIRE_NOTHROW( withName::get("ETHANOL")              );
    REQUIRE_NOTHROW( withName::get("1-PROPANOL")           );
    REQUIRE_NOTHROW( withName::get("1-BUTANOL")            );
    REQUIRE_NOTHROW( withName::get("1-HEXANOL")            );
    REQUIRE_NOTHROW( withName::get("2-PROPANOL")           );
    REQUIRE_NOTHROW( withName::get("PHENOL")               );
    REQUIRE_NOTHROW( withName::get("ETHYLENE-GLYCOL")      );
    REQUIRE_NOTHROW( withName::get("ACETIC-ACID")          );
    REQUIRE_NOTHROW( withName::get("N-BUTYRIC-ACID")       );
    REQUIRE_NOTHROW( withName::get("BENZOIC-ACID")         );
    REQUIRE_NOTHROW( withName::get("ACETONITRILE")         );
    REQUIRE_NOTHROW( withName::get("METHYLAMINE")          );
    REQUIRE_NOTHROW( withName::get("ETHYLAMINE")           );
    REQUIRE_NOTHROW( withName::get("NITROmethane")         );
    REQUIRE_NOTHROW( withName::get("CARBOn-TETRACHLORIDE") );
    REQUIRE_NOTHROW( withName::get("CHLORoform")           );
    REQUIRE_NOTHROW( withName::get("DICHLOROmethane")      );
    REQUIRE_NOTHROW( withName::get("METHYL-Chloride")      );
    REQUIRE_NOTHROW( withName::get("ETHYL-CHloride")       );
    REQUIRE_NOTHROW( withName::get("CHLOROBENZENE")        );
    REQUIRE_NOTHROW( withName::get("TETRAFLUOROETHANE")    );
    REQUIRE_NOTHROW( withName::get("ARGON")                );
    REQUIRE_NOTHROW( withName::get("KRYPTON")              );
    REQUIRE_NOTHROW( withName::get("XENON")                );
    REQUIRE_NOTHROW( withName::get("HELIUM")               );
    REQUIRE_NOTHROW( withName::get("HYDROGEN")             );
    REQUIRE_NOTHROW( withName::get("OXYGEN")               );
    REQUIRE_NOTHROW( withName::get("NITROGEN")             );
    REQUIRE_NOTHROW( withName::get("AIR")                  );
    REQUIRE_NOTHROW( withName::get("CHLORINE")             );
    REQUIRE_NOTHROW( withName::get("CARBON-MONOXIDE")      );
    REQUIRE_NOTHROW( withName::get("CARBON-DIOXIDE")       );
    REQUIRE_NOTHROW( withName::get("CARBON-DISULFIDE")     );
    REQUIRE_NOTHROW( withName::get("HYDROGEN-SULFIDE")     );
    REQUIRE_NOTHROW( withName::get("SULFUR-DIOXIDE")       );
    REQUIRE_NOTHROW( withName::get("SULFUR-TRIOXIDE")      );
    REQUIRE_NOTHROW( withName::get("NITRIC-OXIDE")         );
    REQUIRE_NOTHROW( withName::get("NITROUS-OXIDE")        );
    REQUIRE_NOTHROW( withName::get("HYDROGEN-CHLORIDE")    );
    REQUIRE_NOTHROW( withName::get("HYDROGEN-CYANIDE")     );
    REQUIRE_NOTHROW( withName::get("WATER")                );
    REQUIRE_NOTHROW( withName::get("AMMONIA")              );
    REQUIRE_NOTHROW( withName::get("NITRIC-ACID")          );
    REQUIRE_NOTHROW( withName::get("SULFURIC-ACID")        );
    REQUIRE_NOTHROW( withName::get("RADON")                );
    REQUIRE_NOTHROW( withName::get("NEON")                 );

    // Checking if get is indeed case insensitive
    REQUIRE_NOTHROW( withName::get("methane") );
    REQUIRE_NOTHROW( withName::get("ethane")  );
    REQUIRE_NOTHROW( withName::get("propane") );

    // Checking if get finds species with other names
    REQUIRE_NOTHROW( withName::get("CH4")   );
    REQUIRE_NOTHROW( withName::get("C2H6")  );
    REQUIRE_NOTHROW( withName::get("C3H8")  );
    REQUIRE_NOTHROW( withName::get("C5H12") );
    REQUIRE_NOTHROW( withName::get("C6H14") );
    REQUIRE_NOTHROW( withName::get("C6H12") );
    REQUIRE_NOTHROW( withName::get("C2H4")  );
    REQUIRE_NOTHROW( withName::get("C3H6")  );
    REQUIRE_NOTHROW( withName::get("C4H8")  );
    REQUIRE_NOTHROW( withName::get("C2H2")  );
    REQUIRE_NOTHROW( withName::get("C6H6")  );
    REQUIRE_NOTHROW( withName::get("C7H8")  );
    REQUIRE_NOTHROW( withName::get("CH4O")  );
    REQUIRE_NOTHROW( withName::get("C2H6O") );
    REQUIRE_NOTHROW( withName::get("Ar")    );
    REQUIRE_NOTHROW( withName::get("Kr")    );
    REQUIRE_NOTHROW( withName::get("Xe")    );
    REQUIRE_NOTHROW( withName::get("He")    );
    REQUIRE_NOTHROW( withName::get("H2")    );
    REQUIRE_NOTHROW( withName::get("O2")    );
    REQUIRE_NOTHROW( withName::get("N2")    );
    REQUIRE_NOTHROW( withName::get("Cl2")   );
    REQUIRE_NOTHROW( withName::get("CO")    );
    REQUIRE_NOTHROW( withName::get("CO2")   );
    REQUIRE_NOTHROW( withName::get("CS2")   );
    REQUIRE_NOTHROW( withName::get("H2S")   );
    REQUIRE_NOTHROW( withName::get("SO2")   );
    REQUIRE_NOTHROW( withName::get("SO3")   );
    REQUIRE_NOTHROW( withName::get("NO")    );
    REQUIRE_NOTHROW( withName::get("N2O")   );
    REQUIRE_NOTHROW( withName::get("HCl")   );
    REQUIRE_NOTHROW( withName::get("HCN")   );
    REQUIRE_NOTHROW( withName::get("H2O")   );
    REQUIRE_NOTHROW( withName::get("NH3")   );
    REQUIRE_NOTHROW( withName::get("HNO3")  );
    REQUIRE_NOTHROW( withName::get("H2SO4") );
    REQUIRE_NOTHROW( withName::get("Ra")    );
    REQUIRE_NOTHROW( withName::get("Ne")    );

    REQUIRE_NOTHROW( withNames::get({ "methane", "CH4" }) );
    REQUIRE_NOTHROW( withNames::get({ "argon", "Ar" }) );
    REQUIRE_THROWS( withNames::get({ "Ab", "Cd" }) );

    REQUIRE( withName::temperature("CO2")    == Approx(304.20)   );
    REQUIRE( withName::pressure("CO2")       == Approx(73.83e+5) );
    REQUIRE( withName::acentricFactor("CO2") == Approx(0.2240)   );

    REQUIRE_THROWS( CriticalProps::append({ { 400.0, 100.0e+5, 0.1234 }, {"HCl"} }) ); // trying to append new entry with same name as an existing one

    CriticalProps::overwrite({ { 400.0, 100.0e+5, 0.1234 }, {"HCl"} });

    REQUIRE( withName::temperature("HCl")    == Approx(400.0)    );
    REQUIRE( withName::pressure("HCl")       == Approx(100.0e+5) );
    REQUIRE( withName::acentricFactor("HCl") == Approx(0.1234)   );

    SubstanceCriticalProps crprops({"HCl", "HYDROCHLORIC-ACID"});
    crprops.setTemperature(200, "celsius");
    crprops.setPressure(150, "bar");
    crprops.setAcentricFactor(0.9999);

    CriticalProps::overwrite(crprops); // overwrite instead of append if there is an existing substance with one or more common names

    REQUIRE( withName::temperature("HCl")    == Approx(473.15)   );
    REQUIRE( withName::pressure("HCl")       == Approx(150.0e+5) );
    REQUIRE( withName::acentricFactor("HCl") == Approx(0.9999)   );

    REQUIRE_NOTHROW( CriticalProps::setMissingAs("H2") ); // set default properties for missing substances as that of H2

    REQUIRE( withName::temperature("InexistentSubstance") == withName::temperature("H2") );
    REQUIRE( withName::pressure("InexistentSubstance") == withName::pressure("H2") );
    REQUIRE( withName::acentricFactor("InexistentSubstance") == withName::acentricFactor("H2") );

    REQUIRE_THROWS( CriticalProps::setMissingAs("XYZ") ); // cannot set default properties for missing substances using inexistent substance in the database
}
