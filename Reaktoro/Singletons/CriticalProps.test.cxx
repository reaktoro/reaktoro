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
#include <Reaktoro/Singletons/CriticalProps.hpp>
using namespace Reaktoro;

namespace withFormula {

auto get(std::string formula) { return CriticalProps::getWithFormula(formula).value(); }
auto temperature(std::string formula) { return get(formula).temperature(); }
auto pressure(std::string formula) { return get(formula).pressure(); }
auto acentricFactor(std::string formula) { return get(formula).acentricFactor(); }

} // withFormula

namespace withName {

auto get(std::string name) { return CriticalProps::getWithName(name).value(); }
auto temperature(std::string name) { return get(name).temperature(); }
auto pressure(std::string name) { return get(name).pressure(); }
auto acentricFactor(std::string name) { return get(name).acentricFactor(); }

} // withFormula

TEST_CASE("Testing CriticalProps", "[CriticalProps]")
{
    REQUIRE(CriticalProps::size() == CriticalProps::substances().size());

    REQUIRE_NOTHROW( withName::get("Argon")            );
    REQUIRE_NOTHROW( withName::get("Ethylene")         );
    REQUIRE_NOTHROW( withName::get("Phenol")           );
    REQUIRE_NOTHROW( withName::get("Meta-Cresol")      );
    REQUIRE_NOTHROW( withName::get("Ortho-Cresol")     );
    REQUIRE_NOTHROW( withName::get("Para-Cresol")      );
    REQUIRE_NOTHROW( withName::get("Methane")          );
    REQUIRE_NOTHROW( withName::get("Carbon-Monoxide")  );
    REQUIRE_NOTHROW( withName::get("Carbon-Dioxide")   );
    REQUIRE_NOTHROW( withName::get("Hydrogen")         );
    REQUIRE_NOTHROW( withName::get("Steam")            );
    REQUIRE_NOTHROW( withName::get("Hydrogen-Sulfide") );
    REQUIRE_NOTHROW( withName::get("Helium")           );
    REQUIRE_NOTHROW( withName::get("Krypton")          );
    REQUIRE_NOTHROW( withName::get("Nitrogen")         );
    REQUIRE_NOTHROW( withName::get("Nitrous-Oxide")    );
    REQUIRE_NOTHROW( withName::get("Neon")             );
    REQUIRE_NOTHROW( withName::get("Ammonia")          );
    REQUIRE_NOTHROW( withName::get("Nitric-Oxide")     );
    REQUIRE_NOTHROW( withName::get("Oxygen")           );
    REQUIRE_NOTHROW( withName::get("Radon")            );
    REQUIRE_NOTHROW( withName::get("Sulfur")           );
    REQUIRE_NOTHROW( withName::get("Sulfur-Dioxide")   );
    REQUIRE_NOTHROW( withName::get("Xenon")            );

    REQUIRE_NOTHROW( withFormula::get("Ar")    );
    REQUIRE_NOTHROW( withFormula::get("C2H4")  );
    REQUIRE_NOTHROW( withFormula::get("C6H6O") );
    REQUIRE_NOTHROW( withFormula::get("C7H8O") );
    REQUIRE_NOTHROW( withFormula::get("C7H8O") );
    REQUIRE_NOTHROW( withFormula::get("C7H8O") );
    REQUIRE_NOTHROW( withFormula::get("CH4")   );
    REQUIRE_NOTHROW( withFormula::get("CO")    );
    REQUIRE_NOTHROW( withFormula::get("CO2")   );
    REQUIRE_NOTHROW( withFormula::get("H2")    );
    REQUIRE_NOTHROW( withFormula::get("H2O")   );
    REQUIRE_NOTHROW( withFormula::get("H2S")   );
    REQUIRE_NOTHROW( withFormula::get("He")    );
    REQUIRE_NOTHROW( withFormula::get("Kr")    );
    REQUIRE_NOTHROW( withFormula::get("N2")    );
    REQUIRE_NOTHROW( withFormula::get("N2O")   );
    REQUIRE_NOTHROW( withFormula::get("Ne")    );
    REQUIRE_NOTHROW( withFormula::get("NH3")   );
    REQUIRE_NOTHROW( withFormula::get("NO")    );
    REQUIRE_NOTHROW( withFormula::get("O2")    );
    REQUIRE_NOTHROW( withFormula::get("Rn")    );
    REQUIRE_NOTHROW( withFormula::get("S2")    );
    REQUIRE_NOTHROW( withFormula::get("SO2")   );
    REQUIRE_NOTHROW( withFormula::get("Xe")    );

    REQUIRE( withFormula::temperature("CO2")    == Approx(304.20)   );
    REQUIRE( withFormula::pressure("CO2")       == Approx(73.83e+5) );
    REQUIRE( withFormula::acentricFactor("CO2") == Approx(0.2240)   );

    CriticalProps::append({ "HCl", "HCl", { 400.0, 100.0e+5, 0.1234 } });

    REQUIRE( withFormula::temperature("HCl")    == Approx(400.0)    );
    REQUIRE( withFormula::pressure("HCl")       == Approx(100.0e+5) );
    REQUIRE( withFormula::acentricFactor("HCl") == Approx(0.1234)   );

    withFormula::get("HCl").setTemperature(200, "celsius");
    withFormula::get("HCl").setPressure(150, "bar");
    withFormula::get("HCl").setAcentricFactor(0.9999);

    REQUIRE( withFormula::temperature("HCl")    == Approx(473.15)    );
    REQUIRE( withFormula::pressure("HCl")       == Approx(150.0e+5) );
    REQUIRE( withFormula::acentricFactor("HCl") == Approx(0.9999)   );
}
