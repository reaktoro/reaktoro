// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Extensions/ThermoFun/ThermoFunDatabase.hpp>

using namespace Reaktoro;

namespace test
{

Map<String, Any> databases;

enum DatabaseType
{
    Phreeqc,    // databases defined by the dat files
    Supcrt,     // databases defined by the xml files with names `supcrt**`
    ThermoFun,  // databases defined by the json files with names `**-thermofun`
    Reaktoro    // databases defined by the json or yaml files `**tbl`
};

auto getDatabaseType(const String& string) -> DatabaseType
{
    if(string.find(".dat") < string.size())             // name of the file `.dat`
        return DatabaseType::Phreeqc;
    else if (string.find("supcrt") < string.size())     // name of the file `supcrt`
        return DatabaseType::Supcrt;
    else if (string.find("thermofun") < string.size())  // name of the file `thermofun`
        return DatabaseType::ThermoFun;
    else
        return DatabaseType::Reaktoro;
}

auto initializeDatabases(const String& string) -> void
{
    if(!databases[string].has_value())
    {
        switch (getDatabaseType(string))
        {
            case DatabaseType::Phreeqc:
                databases[string] = PhreeqcDatabase(string);
                break;
            case DatabaseType::Supcrt or DatabaseType::Reaktoro:
                databases[string] = SupcrtDatabase(string);
                break;
            case DatabaseType::ThermoFun:
                databases[string] = ThermoFunDatabase(string);
                break;
            default:
                warning(true, "Reaktoro is not able to understand the type of database in the file " + string);
        }
    }
}
}


//=================================================================================================
// AUXILIARY FUNCTIONS: DECLARATION
//=================================================================================================

/// Return the contents of extra data, in PHREEQC format, as a string, to complement an existing database.
auto getStringContentsPhreeqcDatabaseComplement() -> String;

/// Calculate the actual equilibrium constant of the formation reaction of a product species with given name.
auto lgK(const PhreeqcDatabase& db, real T, real P, String name) -> real;

/// Calculate the actual enthalpy of formation reaction of a product species with given name.
auto dH0(const PhreeqcDatabase& db, real T, real P, String name) -> real;

/// Calculate the actual standard molar Gibbs energies of species with given name.
auto G0(const PhreeqcDatabase& db, real T, real P, String name) -> real;

/// Calculate the expected equilibrium constant correction using van't Hoff equation.
/// @param T The temperature for the calculation (in K)
/// @param log_k0 The equilibrium constant at 298.15K (log base 10)
/// @param delta_h0 The enthalpy of reaction at 298.15K (in J/mol)
auto lgK_vantHoff(real T, real log_k0, real delta_h0) -> real;

/// Calculate the expected equilibrium constant correction using the analytic equation used in PHREEQC.
/// @param T The temperature for the calculation (in K)
/// @param log_k0 The equilibrium constant at 298.15K (log base 10)
/// @param delta_h0 The enthalpy of reaction at 298.15K (in J/mol)
auto lgK_analytic(real T, const Vec<real>& A) -> real;

/// Return -1 for gases and minerals, 1 for other species.
/// This method is needed because Reaktoro considers all reactions in a
/// formation reaction format (i.e., the species associated with the reaction
/// is on the right-hand side; it's a product species). However, PHREEQC
/// considers reactions for gases and minerals in a dissociation format (i.e.,
/// the species associated with the reaction is on the left-hand side; it's a
/// reactant species).
auto lgK_sign(const PhreeqcDatabase& db, String name) -> real;

/// Calculate the expected equilibrium constant of a reaction using standard molar Gibbs energies of species.
auto lgK_fromG0(const PhreeqcDatabase& db, real T, real P, Pairs<String, double> reaction) -> real;

TEST_CASE("Testing PhreeqcDatabase constructor and load method", "[PhreeqcDatabase]")
{
    //-------------------------------------------------------------------------
    // Testing constructor with name of embedded database file
    //-------------------------------------------------------------------------
    CHECK_NOTHROW( PhreeqcDatabase("phreeqc.dat") );
    CHECK_NOTHROW( PhreeqcDatabase("pitzer.dat") );

    CHECK_THROWS( PhreeqcDatabase("xyz.dat") ); // there is no xyz.dat embedded database!

    //-------------------------------------------------------------------------
    // Testing method PhreeqcDatabase::load
    //-------------------------------------------------------------------------
    PhreeqcDatabase db;

    // Load using a string containing the database contents in PHREEQC format
    db.load(PhreeqcDatabase::contents("pitzer.dat"));

    // Checking the existence of elements in the customized database
    CHECK_NOTHROW( db.elements().index("H")  );
    CHECK_NOTHROW( db.elements().index("C")  );
    CHECK_NOTHROW( db.elements().index("O")  );
    CHECK_NOTHROW( db.elements().index("Na") );
    CHECK_NOTHROW( db.elements().index("Cl") );
    CHECK_NOTHROW( db.elements().index("Ca") );

    // Checking the inexistence of element U (which will be added later with a database complement)
    CHECK_THROWS( db.elements().index("U")  );

    // Checking the existence of species in the customized database
    CHECK_NOTHROW( db.species().index("H2O")   );
    CHECK_NOTHROW( db.species().index("H+")    );
    CHECK_NOTHROW( db.species().index("OH-")   );
    CHECK_NOTHROW( db.species().index("CO2")   );
    CHECK_NOTHROW( db.species().index("CO3-2") );
    CHECK_NOTHROW( db.species().index("HCO3-") );

    // Checking the inexistence of species with element U (which will be added later with a database complement)
    CHECK_THROWS( db.species().index("U+4")           );
    CHECK_THROWS( db.species().index("UO2+")          );
    CHECK_THROWS( db.species().index("UO2+2")         );
    CHECK_THROWS( db.species().index("U(OH)4")        );
    CHECK_THROWS( db.species().index("U(OH)5-")       );
    CHECK_THROWS( db.species().index("UO2OH+")        );
    CHECK_THROWS( db.species().index("(UO2)2(OH)2+2") );
    CHECK_THROWS( db.species().index("(UO2)3(OH)5+")  );
    CHECK_THROWS( db.species().index("UO2CO3")        );
    CHECK_THROWS( db.species().index("UO2(CO3)2-2")   );
    CHECK_THROWS( db.species().index("UO2(CO3)3-4")   );

    // // // In this new load operation, complement the PhreeqcDatabase object with more contents
    // db.extend(getStringContentsPhreeqcDatabaseComplement()); // TODO Implement PhreeqcDatabase::extend method using PhreeqcUtis::execute for this.

    // // Checking the existence of newly added element U
    // CHECK_NOTHROW( db.elements().index("U")  );

    // // Checking the existence of newly added species with element U
    // CHECK_NOTHROW( db.species().index("U+4")           );
    // CHECK_NOTHROW( db.species().index("UO2+")          );
    // CHECK_NOTHROW( db.species().index("UO2+2")         );
    // CHECK_NOTHROW( db.species().index("U(OH)4")        );
    // CHECK_NOTHROW( db.species().index("U(OH)5-")       );
    // CHECK_NOTHROW( db.species().index("UO2OH+")        );
    // CHECK_NOTHROW( db.species().index("(UO2)2(OH)2+2") );
    // CHECK_NOTHROW( db.species().index("(UO2)3(OH)5+")  );
    // CHECK_NOTHROW( db.species().index("UO2CO3")        );
    // CHECK_NOTHROW( db.species().index("UO2(CO3)2-2")   );
    // CHECK_NOTHROW( db.species().index("UO2(CO3)3-4")   );

    // // Ensure previously existing elements are still present
    // CHECK_NOTHROW( db.elements().index("H")  );
    // CHECK_NOTHROW( db.elements().index("C")  );
    // CHECK_NOTHROW( db.elements().index("O")  );
    // CHECK_NOTHROW( db.elements().index("Na") );
    // CHECK_NOTHROW( db.elements().index("Cl") );
    // CHECK_NOTHROW( db.elements().index("Ca") );

    // // Ensure previously existing species are still present
    // CHECK_NOTHROW( db.species().index("H2O")   );
    // CHECK_NOTHROW( db.species().index("H+")    );
    // CHECK_NOTHROW( db.species().index("OH-")   );
    // CHECK_NOTHROW( db.species().index("CO2")   );
    // CHECK_NOTHROW( db.species().index("CO3-2") );
    // CHECK_NOTHROW( db.species().index("HCO3-") );
}

TEST_CASE("Testing species and its attributes after constructing PhreeqcDatabase", "[PhreeqcDatabase]")
{
    const PhreeqcDatabase db("phreeqc.dat");

    //-------------------------------------------------------------------------
    // Testing for the existence of selected elements
    //-------------------------------------------------------------------------

    CHECK_NOTHROW( db.elements().index("H")  );
    CHECK_NOTHROW( db.elements().index("C")  );
    CHECK_NOTHROW( db.elements().index("O")  );
    CHECK_NOTHROW( db.elements().index("Na") );
    CHECK_NOTHROW( db.elements().index("Cl") );
    CHECK_NOTHROW( db.elements().index("K")  );
    CHECK_NOTHROW( db.elements().index("Ca") );
    CHECK_NOTHROW( db.elements().index("Mg") );
    CHECK_NOTHROW( db.elements().index("Si") );
    CHECK_NOTHROW( db.elements().index("Fe") );

    //-------------------------------------------------------------------------
    // Testing for the existence of selected species
    //-------------------------------------------------------------------------
    CHECK_NOTHROW( db.species().index("H2O")         );
    CHECK_NOTHROW( db.species().index("H+")          );
    CHECK_NOTHROW( db.species().index("CO3-2")       );
    CHECK_NOTHROW( db.species().index("Calcite")     );
    CHECK_NOTHROW( db.species().index("Quartz")      );
    CHECK_NOTHROW( db.species().index("CO2(g)")      );
    CHECK_NOTHROW( db.species().index("NaX")         );
    CHECK_NOTHROW( db.species().index("X-")          );
    CHECK_NOTHROW( db.species().index("Hfo_sOHCa+2") );

    Species species;

    species = db.species().get("H2O");

    CHECK( species.name() == "H2O"                               );
    CHECK( species.formula().str() == "H2O"                      );
    CHECK( species.substance() == "H2O"                          );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(0.0180160)              );
    CHECK( species.elements().size() == 2                        );
    CHECK( species.elements().coefficient("H") == 2              );
    CHECK( species.elements().coefficient("O") == 1              );
    CHECK( species.reaction().reactants().size() == 0            );

    species = db.species().get("e-");

    CHECK( species.name() == "e-"                                );
    CHECK( species.formula().str() == "e-"                       );
    CHECK( species.substance() == "e-"                           );
    CHECK( species.charge() == -1                                );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(5.4857990888E-07)       );
    CHECK( species.elements().size() == 1                        );
    CHECK( species.elements().coefficient("e") == 1              );
    CHECK( species.reaction().reactants().size() == 0            );

    species = db.species().get("CO2");

    CHECK( species.name() == "CO2"                               );
    CHECK( species.formula().str() == "CO2"                      );
    CHECK( species.substance() == "CO2"                          );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(0.0440111)              );
    CHECK( species.elements().size() == 2                        );
    CHECK( species.elements().coefficient("C") == 1              );
    CHECK( species.elements().coefficient("O") == 2              );
    CHECK( species.reaction().reactants().size() == 3            );
    CHECK( species.reaction().stoichiometry("CO3-2") == 1        );
    CHECK( species.reaction().stoichiometry("H+") == 2           );
    CHECK( species.reaction().stoichiometry("H2O") == -1         );

    species = db.species().get("CO3-2");

    CHECK( species.name() == "CO3-2"                             );
    CHECK( species.formula().str() == "CO3-2"                    );
    CHECK( species.substance() == "CO3-2"                        );
    CHECK( species.charge() == -2                                );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(0.0600121972)           );
    CHECK( species.elements().size() == 2                        );
    CHECK( species.elements().coefficient("C") == 1              );
    CHECK( species.elements().coefficient("O") == 3              );
    CHECK( species.reaction().reactants().size() == 0            );

    species = db.species().get("Pb2OH+3");

    CHECK( species.name() == "Pb2OH+3"                           );
    CHECK( species.formula().str() == "Pb2OH+3"                  );
    CHECK( species.substance() == "Pb2OH+3"                      );
    CHECK( species.charge() == +3                                );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(0.4313880)              );
    CHECK( species.elements().size() == 3                        );
    CHECK( species.elements().coefficient("Pb") == 2             );
    CHECK( species.elements().coefficient("O") == 1              );
    CHECK( species.elements().coefficient("H") == 1              );
    CHECK( species.reaction().reactants().size() == 3            );
    CHECK( species.reaction().stoichiometry("Pb+2") == 2         );
    CHECK( species.reaction().stoichiometry("H2O") == 1          );
    CHECK( species.reaction().stoichiometry("H+") == -1          );

    species = db.species().get("CO2(g)");

    CHECK( species.name() == "CO2(g)"                            );
    CHECK( species.formula().str() == "CO2"                      );
    CHECK( species.substance() == "CO2"                          );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Gas       );
    CHECK( species.molarMass() == Approx(0.0440111)              );
    CHECK( species.elements().size() == 2                        );
    CHECK( species.elements().coefficient("C") == 1              );
    CHECK( species.elements().coefficient("O") == 2              );
    CHECK( species.reaction().reactants().size() == 1            );
    CHECK( species.reaction().stoichiometry("CO2") == 1          );

    species = db.species().get("H2S(g)");

    CHECK( species.name() == "H2S(g)"                            );
    CHECK( species.formula().str() == "H2S"                      );
    CHECK( species.substance() == "H2S"                          );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Gas       );
    CHECK( species.molarMass() == Approx(0.0340800)              );
    CHECK( species.elements().size() == 2                        );
    CHECK( species.elements().coefficient("H") == 2              );
    CHECK( species.elements().coefficient("S") == 1              );
    CHECK( species.reaction().reactants().size() == 2            );
    CHECK( species.reaction().stoichiometry("H+") == 1           );
    CHECK( species.reaction().stoichiometry("HS-") == 1          );

    species = db.species().get("Calcite");

    CHECK( species.name() == "Calcite"                           );
    CHECK( species.formula().str() == "CaCO3"                    );
    CHECK( species.substance() == "Calcite"                      );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Solid     );
    CHECK( species.molarMass() == Approx(0.1000911)              );
    CHECK( species.elements().size() == 3                        );
    CHECK( species.elements().coefficient("Ca") == 1             );
    CHECK( species.elements().coefficient("C") == 1              );
    CHECK( species.elements().coefficient("O") == 3              );
    CHECK( species.reaction().reactants().size() == 2            );
    CHECK( species.reaction().stoichiometry("Ca+2") == 1         );
    CHECK( species.reaction().stoichiometry("CO3-2") == 1        );

    species = db.species().get("Hydroxyapatite");

    CHECK( species.name() == "Hydroxyapatite"                    );
    CHECK( species.formula().str() == "Ca5(PO4)3OH"              );
    CHECK( species.substance() == "Hydroxyapatite"               );
    CHECK( species.charge() == 0                                 );
    CHECK( species.aggregateState() == AggregateState::Solid     );
    CHECK( species.molarMass() == Approx(0.5023294)              );
    CHECK( species.elements().size() == 4                        );
    CHECK( species.elements().coefficient("Ca") == 5             );
    CHECK( species.elements().coefficient("P") == 3              );
    CHECK( species.elements().coefficient("O") == 13             );
    CHECK( species.elements().coefficient("H") == 1              );
    CHECK( species.reaction().reactants().size() == 4            );
    CHECK( species.reaction().stoichiometry("H2O") == 1          );
    CHECK( species.reaction().stoichiometry("HPO4-2") == 3       );
    CHECK( species.reaction().stoichiometry("Ca+2") == 5         );
    CHECK( species.reaction().stoichiometry("H+") == -4          );
}

TEST_CASE("Testing standard thermodynamic properties calculations", "[PhreeqcDatabase]")
{
    static const PhreeqcDatabase db("phreeqc.dat"); // if not static, GENERATE causes many instantiations of PhreeqcDatabase

    // const auto T = GENERATE( 298.15, 303.15, 333.15, 363.15 );
    // const auto P = GENERATE( 1.0e5, 10.0e5, 50.0e5, 100.0e5 );
    const auto T = GENERATE( 298.15 );
    const auto P = GENERATE( 1.0e5 );

    //---------------------------------------------------------------------------------------------
    // Testing temperature correction of equilibrium constants using van't Hoff equation
    //---------------------------------------------------------------------------------------------

    // The conversion constant from kcal to J.
    const auto kcal_to_joule = 4184.0;

    // The equilibrium constant (log_k) and enthalpy of reaction (delta_h) (in J) of some selected species
    const Map<String, Pair<real, real>> vanthoff_data =
    {
        { "SrSO4"    , {  2.29  ,  2.08   * kcal_to_joule} },
        { "Cu+"      , {  2.72  ,  1.65   * kcal_to_joule} },
        { "CuCl2-"   , {  5.50  , -0.42   * kcal_to_joule} },
        { "CuCO3"    , {  6.73  ,  0.00   * kcal_to_joule} },
        { "Zn(OH)4-2", { -41.2  ,  0.00   * kcal_to_joule} },
        { "Kaolinite", {  7.435 , -35.300 * kcal_to_joule} },
        { "Albite"   , { -18.002,  25.896 * kcal_to_joule} },
    };

    for(auto [species, pair] : vanthoff_data)
    {
        INFO("species: " << species);
        INFO("T: " << T << " K");
        INFO("P: " << P << " Pa");
        CHECK( lgK(db, T, P, species) == Approx(lgK_vantHoff(T, pair.first, pair.second) * lgK_sign(db, species)) );
    }

    //---------------------------------------------------------------------------------------------
    // Testing temperature correction of equilibrium constants using analytic expression
    //---------------------------------------------------------------------------------------------

    // The analytic coefficients for temperature correction of logK of some selected species
    const Map<String, Vec<real>> coefficients_data =
    {   // Name         A1, A2, A3, A4, A5, A6
        { "CO2"    , {{ 464.1965, 0.09344813, -26986.16, -165.75951, 2248628.9, 0 }} },
        { "CaCO3"  , {{ -1228.732, -0.299440, 35512.75, 485.818, 0, 0 }} },
        { "SrCO3"  , {{ -1.019, 0.012826, 0, 0, 0, 0 }} },
        { "Calcite", {{ -171.9065, -0.077993, 2839.319, 71.595, 0, 0 }} },
    };

    for(auto [species, A] : coefficients_data)
    {
        INFO("species: " << species);
        INFO("T: " << T << " K");
        INFO("P: " << P << " Pa");
        CHECK( lgK(db, T, P, species) == Approx(lgK_analytic(T, A) * lgK_sign(db, species)) );
    }

    // The analytic coefficients for temperature correction of logK of some selected species
    const Map<String, Pairs<String, double>> reactions =
    {   // Species   Reaction
        { "OH-"    , parseReactionEquation("H2O = OH- + H+")           },
        { "CO2"    , parseReactionEquation("CO3-2 + 2*H+ = CO2 + H2O") },
        { "CaCO3"  , parseReactionEquation("Ca+2 + CO3-2 = CaCO3")     },
        { "SrCO3"  , parseReactionEquation("Sr+2 + CO3-2 = SrCO3")     },
        { "CO2(g)" , parseReactionEquation("CO2 = CO2(g)")             },
        { "Calcite", parseReactionEquation("CO3-2 + Ca+2 = Calcite")   },
    };

    for(auto [species, reaction] : reactions)
    {
        INFO("species: " << species);
        INFO("T: " << T << " K");
        INFO("P: " << P << " Pa");
        CHECK( lgK(db, T, P, species) == Approx(lgK_fromG0(db, T, P, reaction)) );
    }
}

TEST_CASE("Testing pressure correction in standard thermodynamic properties calculations", "[PhreeqcDatabase]")
{
    test::initializeDatabases("phreeqc.dat");
    static const auto db = std::any_cast<PhreeqcDatabase>(test::databases["phreeqc.dat"]);

    const auto T = GENERATE( 298.15 );
    const auto P = GENERATE( 1.0e5 );

    const auto T0 = 298.15;   // in K (equivalent to 25 C)
    const auto P0 = 101325.0; // in Pa (equivalent to 1 atm)

    const auto T1 = 298.15 + 60.0; // in K (equivalent to 85 C)
    const auto P1 = 101325.0 * 10; // in Pa (equivalent to 10 atm)

    {
        const auto species = db.species().get("H2O");

        const auto propsT0P0 = species.props(T0, P0);
        const auto propsT1P0 = species.props(T1, P0);
        const auto propsT0P1 = species.props(T0, P1);
        const auto propsT1P1 = species.props(T1, P1);

        INFO("species: " << species.name());

        CHECK(propsT0P0.G0  == Approx(0.0)         );
        CHECK(propsT0P0.H0  == Approx(0.0)         );
        CHECK(propsT0P0.V0  == Approx(1.80694e-05) );
        CHECK(propsT0P0.Cp0 == Approx(0.0)         );

        CHECK(propsT1P0.G0  == Approx(0.0)         );
        CHECK(propsT1P0.H0  == Approx(0.0)         );
        CHECK(propsT1P0.V0  == Approx(1.85999e-05) );
        CHECK(propsT1P0.Cp0 == Approx(0.0)         );

        CHECK(propsT0P1.G0  == Approx(0.0)         );
        CHECK(propsT0P1.H0  == Approx(0.0)         );
        CHECK(propsT0P1.V0  == Approx(1.80621e-05) );
        CHECK(propsT0P1.Cp0 == Approx(0.0)         );

        CHECK(propsT1P1.G0  == Approx(0.0)         );
        CHECK(propsT1P1.H0  == Approx(0.0)         );
        CHECK(propsT1P1.V0  == Approx(1.85919e-05) );
        CHECK(propsT1P1.Cp0 == Approx(0.0)         );
    }

    {
        const auto species = db.species().get("H+");

        const auto propsT0P0 = species.props(T0, P0);
        const auto propsT1P0 = species.props(T1, P0);
        const auto propsT0P1 = species.props(T0, P1);
        const auto propsT1P1 = species.props(T1, P1);

        INFO("species: " << species.name());

        CHECK(propsT0P0.G0  == Approx(0.0) );
        CHECK(propsT0P0.H0  == Approx(0.0) );
        CHECK(propsT0P0.V0  == Approx(0.0) );
        CHECK(propsT0P0.Cp0 == Approx(0.0) );

        CHECK(propsT1P0.G0  == Approx(0.0) );
        CHECK(propsT1P0.H0  == Approx(0.0) );
        CHECK(propsT1P0.V0  == Approx(0.0) );
        CHECK(propsT1P0.Cp0 == Approx(0.0) );

        CHECK(propsT0P1.G0  == Approx(0.0) );
        CHECK(propsT0P1.H0  == Approx(0.0) );
        CHECK(propsT0P1.V0  == Approx(0.0) );
        CHECK(propsT0P1.Cp0 == Approx(0.0) );

        CHECK(propsT1P1.G0  == Approx(0.0) );
        CHECK(propsT1P1.H0  == Approx(0.0) );
        CHECK(propsT1P1.V0  == Approx(0.0) );
        CHECK(propsT1P1.Cp0 == Approx(0.0) );
    }

    {
        const auto species = db.species().get("CO2");

        const auto propsT0P0 = species.props(T0, P0);
        const auto propsT1P0 = species.props(T1, P0);
        const auto propsT0P1 = species.props(T0, P1);
        const auto propsT1P1 = species.props(T1, P1);

        INFO("species: " << species.name());

        CHECK(propsT0P0.G0  == Approx(-95213.7)    );
        CHECK(propsT0P0.H0  == Approx(-24010.3)    );
        CHECK(propsT0P0.V0  == Approx(3.44329e-05) );
        CHECK(propsT0P0.Cp0 == Approx(0.0)         );

        CHECK(propsT1P0.G0  == Approx(-113039)     );
        CHECK(propsT1P0.H0  == Approx(12124.6)     );
        CHECK(propsT1P0.V0  == Approx(3.77425e-05) );
        CHECK(propsT1P0.Cp0 == Approx(0.0)         );

        CHECK(propsT0P1.G0  == Approx(-95161)    );
        CHECK(propsT0P1.H0  == Approx(-23957.6)    );
        CHECK(propsT0P1.V0  == Approx(3.44286e-05) );
        CHECK(propsT0P1.Cp0 == Approx(0.0)         );

        CHECK(propsT1P1.G0  == Approx(-112978)     );
        CHECK(propsT1P1.H0  == Approx(12185.4)     );
        CHECK(propsT1P1.V0  == Approx(3.77222e-05) );
        CHECK(propsT1P1.Cp0 == Approx(0.0)         );
    }

    {
        const auto species = db.species().get("Calcite");

        const auto propsT0P0 = species.props(T0, P0);
        const auto propsT1P0 = species.props(T1, P0);
        const auto propsT0P1 = species.props(T0, P1);
        const auto propsT1P1 = species.props(T1, P1);

        INFO("species: " << species.name());

        CHECK(propsT0P0.G0  == Approx(-48402.9) );
        CHECK(propsT0P0.H0  == Approx(9608.99)  );
        CHECK(propsT0P0.V0  == Approx(3.69e-05) );
        CHECK(propsT0P0.Cp0 == Approx(0.0)      );

        CHECK(propsT1P0.G0  == Approx(-62078.3) );
        CHECK(propsT1P0.H0  == Approx(32690.1)  );
        CHECK(propsT1P0.V0  == Approx(3.69e-05) );
        CHECK(propsT1P0.Cp0 == Approx(0.0)      );

        CHECK(propsT0P1.G0  == Approx(-48347.9) );
        CHECK(propsT0P1.H0  == Approx(9664.07)  );
        CHECK(propsT0P1.V0  == Approx(3.69e-05) );
        CHECK(propsT0P1.Cp0 == Approx(0.0)      );

        CHECK(propsT1P1.G0  == Approx(-62017.7) );
        CHECK(propsT1P1.H0  == Approx(32750.7)  );
        CHECK(propsT1P1.V0  == Approx(3.69e-05) );
        CHECK(propsT1P1.Cp0 == Approx(0.0)      );
    }
}

//=================================================================================================
// AUXILIARY FUNCTIONS: DEFINITION
//=================================================================================================

auto evalReactionThermoProps(const Species& species, real T, real P)
{
    const auto dV0 = 0.0;
    return species.reaction().reactionThermoModel()({T, P, dV0});
}

auto lgK(const PhreeqcDatabase& db, real T, real P, String name) -> real
{
    const auto species = db.species().get(name);
    const auto dG0 = evalReactionThermoProps(species, T, P).dG0;
    const auto R = universalGasConstant;
    const auto lnK = -dG0/(R*T);
    return lnK / ln10;
}

auto dH0(const PhreeqcDatabase& db, real T, real P, String name) -> real
{
    const auto species = db.species().get(name);
    const auto dH0 = evalReactionThermoProps(species, T, P).dH0;
    return dH0;
}

auto G0(const PhreeqcDatabase& db, real T, real P, String name) -> real
{
    const auto species = db.species().get(name);
    return species.standardThermoProps(T, P).G0;
}

auto lgK_vantHoff(real T, real log_k0, real delta_h0) -> real
{
    const auto T0 = 298.15;
    const auto R  = universalGasConstant;
    return log_k0 - (delta_h0 * (1.0/T - 1.0/T0)/R) / ln10;
}

auto lgK_analytic(real T, const Vec<real>& A) -> real
{
    assert(A.size() == 6);
    return A[0] + A[1]*T + A[2]/T + A[3]*log10(T) + A[4]/(T*T) + A[5]*T*T;
}

auto lgK_sign(const PhreeqcDatabase& db, String name) -> real
{
    // NOTE: The logic below is needed because in Reaktoro any reaction
    // associated with a chemical species is represented as a formation
    // reaction, i.e., a reaction in which the species it self is on the
    // product side (right-hand side of the reaction). In PHREEQC, this does
    // not for gases and minerals. That's why we multiply the log(K) by -1 for
    // chemical species with gas and solid aggregate state.
    const auto species = db.species().get(name);
    const auto is_phreeqc_phase = oneof(species.aggregateState(), AggregateState::Gas, AggregateState::Solid);
    const auto sign = is_phreeqc_phase ? -1.0 : 1.0;
    return sign;
}

auto lgK_fromG0(const PhreeqcDatabase& db, real T, real P, Pairs<String, double> reaction) -> real
{
    real lnK = 0.0;
    for(auto [species, coeff] : reaction)
        lnK += coeff * G0(db, T, P, species);
    lnK /= -universalGasConstant * T;
    const auto lgK = lnK/ln10;
    return lgK;
}

auto getStringContentsPhreeqcDatabaseComplement() -> String
{
    return R"xyz(
SOLUTION_MASTER_SPECIES
        U       U+4     0.0     238.0290     238.0290
        U(4)    U+4     0.0     238.0290
        U(5)    UO2+    0.0     238.0290
        U(6)    UO2+2   0.0     238.0290
SOLUTION_SPECIES
        #primary master species for U
        #is also secondary master species for U(4)
        U+4 = U+4
                log_k          0.0
        U+4 + 4 H2O = U(OH)4 + 4 H+
                log_k          -8.538
                delta_h        24.760 kcal
        U+4 + 5 H2O = U(OH)5- + 5 H+
                log_k          -13.147
                delta_h        27.580 kcal
        #secondary master species for U(5)
        U+4 + 2 H2O = UO2+ + 4 H+ + e-
                log_k          -6.432
                delta_h        31.130 kcal
        #secondary master species for U(6)
        U+4 + 2 H2O = UO2+2 + 4 H+ + 2 e-
                log_k          -9.217
                delta_h        34.430 kcal
        UO2+2 + H2O = UO2OH+ + H+
                log_k          -5.782
                delta_h        11.015 kcal
        2UO2+2 + 2H2O = (UO2)2(OH)2+2 + 2H+
                log_k          -5.626
                delta_h        -36.04 kcal
        3UO2+2 + 5H2O = (UO2)3(OH)5+ + 5H+
                log_k          -15.641
                delta_h        -44.27 kcal
        UO2+2 + CO3-2 = UO2CO3
                log_k          10.064
                delta_h        0.84 kcal
        UO2+2 + 2CO3-2 = UO2(CO3)2-2
                log_k          16.977
                delta_h        3.48 kcal
        UO2+2 + 3CO3-2 = UO2(CO3)3-4
                log_k          21.397
                delta_h        -8.78 kcal
PHASES
        Uraninite
        UO2 + 4 H+ = U+4 + 2 H2O
        log_k          -3.490
        delta_h        -18.630 kcal
END
)xyz";
}
