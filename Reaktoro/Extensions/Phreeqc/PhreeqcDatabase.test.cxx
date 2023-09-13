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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>

using namespace Reaktoro;

namespace test {

auto getPhreeqcDatabase(const String& name) -> PhreeqcDatabase
{
    // Define map storing names files with corresponding initialized PhreeqcDatabase instances
    static Map<String, PhreeqcDatabase> databases;

    // Check if databases already contains PhreeqcDatabase database corresponding to the input filename
    const auto it = databases.find(name);
    if(it != databases.end()) return it->second;

    // Create PhreeqcDatabase database and add an try to the `databases` map
    PhreeqcDatabase db(name);
    databases[name] = db;

    return db;
}

} // namespace test

//=================================================================================================
// AUXILIARY FUNCTIONS: DECLARATION
//=================================================================================================

/// Return the contents of extra data, in PHREEQC format, as a string, to complement an existing database.
auto getStringContentsPhreeqcDatabaseComplement() -> String;

/// Calculate the actual equilibrium constant of the formation reaction of a product species with given name.
auto lgK(const PhreeqcDatabase& db, real T, real P, const String& name) -> real;

/// Calculate the actual enthalpy of formation reaction of a product species with given name.
auto dH0(const PhreeqcDatabase& db, real T, real P, const String& name) -> real;

/// Calculate the actual isobaric heat capacity of formation reaction of a product species with given name.
auto dCp0(const PhreeqcDatabase& db, real T, real P, const String& name) -> real;

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
auto lgK_analytic(real T, Vec<real> A) -> real;

/// Return -1 for gases and minerals, 1 for other species.
/// This method is needed because Reaktoro considers all reactions in a
/// formation reaction format (i.e., the species associated with the reaction
/// is on the right-hand side; it's a product species). However, PHREEQC
/// considers reactions for gases and minerals in a dissociation format (i.e.,
/// the species associated with the reaction is on the left-hand side; it's a
/// reactant species).
auto lgK_sign(const PhreeqcDatabase& db, const String& name) -> real;

/// Calculate the expected equilibrium constant of a reaction using standard molar Gibbs energies of species.
auto lgK_fromG0(const PhreeqcDatabase& db, real T, real P, const Pairs<String, double>& reaction) -> real;

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

    //-------------------------------------------------------------------------
    // Testing species with same name but different aggregate states are both added
    //-------------------------------------------------------------------------
    {
        PhreeqcDatabase db("llnl.dat");

        auto const aqueous = db.species().withAggregateState(AggregateState::Aqueous);
        auto const solids = db.species().withAggregateState(AggregateState::Solid);

        CHECK_NOTHROW( aqueous.getWithName("Cd(OH)2") );
        CHECK_NOTHROW( aqueous.getWithName("Fe(OH)2") );
        CHECK_NOTHROW( aqueous.getWithName("Fe(OH)3") );
        CHECK_NOTHROW( aqueous.getWithName("FeSO4") );

        CHECK_NOTHROW( solids.getWithName("Cd(OH)2") );
        CHECK_NOTHROW( solids.getWithName("Fe(OH)2") );
        CHECK_NOTHROW( solids.getWithName("Fe(OH)3") );
        CHECK_NOTHROW( solids.getWithName("FeSO4") );
    }

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
    const auto db = test::getPhreeqcDatabase("phreeqc.dat");

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
    CHECK( species.reaction().reactants().empty()                );

    species = db.species().get("e-");

    CHECK( species.name() == "e-"                                );
    CHECK( species.formula().str() == "e-"                       );
    CHECK( species.substance() == "e-"                           );
    CHECK( species.charge() == -1                                );
    CHECK( species.aggregateState() == AggregateState::Aqueous   );
    CHECK( species.molarMass() == Approx(5.4857990888E-07)       );
    CHECK( species.elements().size() == 1                        );
    CHECK( species.elements().coefficient("e") == 1              );
    CHECK( species.reaction().reactants().empty()                );

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
    CHECK( species.reaction().reactants().empty()                );

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
    static const auto db = test::getPhreeqcDatabase("phreeqc.dat");

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
        INFO("species: " << species)
        INFO("T: " << T << " K")
        INFO("P: " << P << " Pa")
        CHECK( lgK(db, T, P, species) == Approx(lgK_vantHoff(T, pair.first, pair.second) * lgK_sign(db, species)) );
    }

    //---------------------------------------------------------------------------------------------
    // Testing temperature correction of equilibrium constants using analytic expression
    //---------------------------------------------------------------------------------------------

    // The analytic coefficients for temperature correction of logK of some selected species
    // IMPORTANT: If any of the tests below start failing, this could be due to
    // an update on phreeqc.dat, which changed the analytical coefficients. In
    // this case, some coefficients below will need to be updated too.
    const Map<String, Vec<real>> coefficients_data =
    {   // Name         A1, A2, A3, A4, A5, A6
        { "OH-" ,         {{ 293.29227, 0.1360833, -10576.913, -123.73158, 0, -6.996455e-5 }} },
        { "HCO3-" ,       {{ 107.8871, 0.03252849, -5151.79, -38.92561, 563713.9 }} },
        { "CO2" ,         {{ 464.1965, 0.09344813, -26986.16, -165.75951, 2248628.9 }} },
        { "(CO2)2" ,      {{ 8.68, -0.0103, -2190 }} },
        { "HSO4-" ,       {{ -56.889, 0.006473, 2307.9, 19.8858 }} },
        { "H2S" ,         {{ -11.17, 0.02386, 3279.0 }} },
        { "HSg-" ,        {{ 11.17, -0.02386, -3279.0 }} },
        { "NH3" ,         {{ 0.6322, -0.001225, -2835.76 }} },
        { "HF" ,          {{ -2.033, 0.012645, 429.01 }} },
        { "CaCO3" ,       {{ -1228.732, -0.299440, 35512.75, 485.818 }} },
        { "CaHCO3+" ,     {{ 1317.0071, 0.34546894, -39916.84, -517.70761, 563713.9 }} },
        { "MgCO3" ,       {{ 0.9910, 0.00667 }} },
        { "MgHCO3+" ,     {{ 48.6721, 0.03252849, -2614.335, -18.00263, 563713.9 }} },
        { "KSO4-" ,       {{ 3.106, 0.0, -673.6 }} },
        { "AlOH+2" ,      {{ -38.253, 0.0, -656.27, 14.327 }} },
        { "Al(OH)2+" ,    {{ 88.50, 0.0, -9391.6, -27.121 }} },
        { "Al(OH)3" ,     {{ 226.374, 0.0, -18247.8, -73.597 }} },
        { "Al(OH)4-" ,    {{ 51.578, 0.0, -11168.9, -14.865 }} },
        { "H3SiO4-" ,     {{ -302.3724, -0.050698, 15669.69, 108.18466, -1119669.0 }} },
        { "H2SiO4-2" ,    {{ -294.0184, -0.072650, 11204.49, 108.18466, -1119669.0 }} },
        { "BaCO3" ,       {{ 0.113, 0.008721 }} },
        { "BaHCO3+" ,     {{ -3.0938, 0.013669 }} },
        { "SrHCO3+" ,     {{ 104.6391, 0.04739549, -5151.79, -38.92561, 563713.9 }} },
        { "SrCO3" ,       {{ -1.019, 0.012826 }} },
        { "Cu2(OH)2+2" ,  {{ 2.497, 0.0, -3833.0 }} },
        // { "Calcite",      {{ -171.9065, -0.077993, 2839.319, 71.595 }} },
        { "Calcite",      {{ 17.118, -0.046528, -3496 }} },
        { "Aragonite",    {{ -171.9773, -0.077993, 2903.293, 71.595 }} },
        { "Strontianite", {{ 155.0305, 0.0, -7239.594, -56.58638 }} },
        { "Witherite",    {{ 607.642, 0.121098, -20011.25, -236.4948 }} },
        { "Gypsum",       {{ 93.7, 5.99E-03, -4e3, -35.019 }} },
        { "Anhydrite",    {{ 84.90, 0, -3135.12, -31.79 }} },
        { "Celestite",    {{ -7.14, 6.11e-3, 75, 0, 0, -1.79e-5 }} },
        { "Barite",       {{ -282.43, -8.972e-2, 5822, 113.08 }} },
        { "Fluorite",     {{ 66.348, 0.0, -4298.2, -25.271 }} },
        { "SiO2(a)",      {{ -0.26, 0.0, -731.0 }} },
        { "Chalcedony",   {{ -0.09, 0.0, -1032.0 }} },
        { "Quartz",       {{ 0.41, 0.0, -1309.0 }} },
        { "Chrysotile",   {{ 13.248, 0.0, 10217.1, -6.1894 }} },
        { "CO2(g)",       {{ 10.5624, -2.3547e-2, -3972.8, 0, 5.8746e5, 1.9194e-5 }} },
        { "H2O(g)",       {{ -16.5066, -2.0013E-3, 2710.7, 3.7646, 0, 2.24E-6 }} },
        { "O2(g)",        {{ -7.5001, 7.8981e-3, 0.0, 0.0, 2.0027e5 }} },
        { "H2(g)",        {{ -9.3114, 4.6473e-3, -49.335, 1.4341, 1.2815e5 }} },
        { "N2(g)",        {{ -58.453, 1.818e-3, 3199, 17.909, -27460 }} },
        // { "H2S(g)",       {{ -97.354, -3.1576e-2, 1.8285e3, 37.44, 28.56 }} },
        { "H2S(g)",       {{ -45.07, -0.02418, 0.0, 17.9205 }} },
        { "CH4(g)",       {{ 10.44, -7.65e-3, -6669, 0, 1.014e6 }} },
        { "NH3(g)",       {{ -18.758, 3.3670e-4, 2.5113e3, 4.8619, 39.192 }} },
        { "Oxg(g)",       {{ -7.5001, 7.8981e-3, 0.0, 0.0, 2.0027e5 }} },
        { "Hdg(g)",       {{ -9.3114, 4.6473e-3, -49.335, 1.4341, 1.2815e5 }} },
        { "Ntg(g)",       {{ -58.453, 1.81800e-3, 3199, 17.909, -27460 }} },
        { "Mtg(g)",       {{ 10.44, -7.65e-3, -6669, 0, 1.014e6 }} },
        { "H2Sg(g)",      {{ -45.07, -0.02418, 0.0, 17.9205 }} },
        { "Melanterite",  {{ 1.447, -0.004153, 0.0, 0.0, -214949.0 }} },
    };

    for(auto [species, A] : coefficients_data)
    {
        INFO("species: " << species)
        INFO("T: " << T << " K")
        INFO("P: " << P << " Pa")
        auto aux1 = lgK(db, T, P, species);
        auto aux2 = lgK_analytic(T, A) * lgK_sign(db, species);
        CHECK( lgK(db, T, P, species) == Approx(lgK_analytic(T, A) * lgK_sign(db, species)) );
    }

    // The analytic coefficients for temperature correction of logK of some selected species
    const Map<String, String> reactions =
    {   // Species   Reaction
        { "OH-"         , "H2O = OH- + H+"                              },
        { "HCO3-"       , "CO3-2 + H+ = HCO3-"                          },
        { "CO2"         , "CO3-2 + 2*H+ = CO2 + H2O"                    },
        { "(CO2)2"      , "2*CO2 = (CO2)2"                              },
        { "HSO4-"       , "SO4-2 + H+ = HSO4-"                          },
        { "H2S"         , "HS- + H+ = H2S"                              },
        { "HSg-"        , "H2Sg = HSg- + H+"                            },
        { "NH3"         , "NH4+ = NH3 + H+"                             },
        { "HF"          , "H+ + F- = HF"                                },
        { "CaCO3"       , "Ca+2 + CO3-2 = CaCO3"                        },
        { "CaHCO3+"     , "Ca+2 + CO3-2 + H+ = CaHCO3+"                 },
        { "MgCO3"       , "Mg+2 + CO3-2 = MgCO3"                        },
        { "MgHCO3+"     , "Mg+2 + H+ + CO3-2 = MgHCO3+"                 },
        { "KSO4-"       , "K+ + SO4-2 = KSO4-"                          },
        { "AlOH+2"      , "Al+3 + H2O = AlOH+2 + H+"                    },
        { "Al(OH)2+"    , "Al+3 + 2*H2O = Al(OH)2+ + 2*H+"              },
        { "Al(OH)3"     , "Al+3 + 3*H2O = Al(OH)3 + 3*H+"               },
        { "Al(OH)4-"    , "Al+3 + 4*H2O = Al(OH)4- + 4*H+"              },
        { "H3SiO4-"     , "H4SiO4 = H3SiO4- + H+"                       },
        { "H2SiO4-2"    , "H4SiO4 = H2SiO4-2 + 2*H+"                    },
        { "BaCO3"       , "Ba+2 + CO3-2 = BaCO3"                        },
        { "BaHCO3+"     , "Ba+2 + HCO3- = BaHCO3+"                      },
        { "SrHCO3+"     , "Sr+2 + CO3-2 + H+ = SrHCO3+"                 },
        { "SrCO3"       , "Sr+2 + CO3-2 = SrCO3"                        },
        { "Cu2(OH)2+2"  , "2*Cu+2 + 2*H2O = Cu2(OH)2+2 + 2*H+"          },
        { "Calcite"     , "Calcite = CO3-2 + Ca+2"                      },
        { "Aragonite"   , "Aragonite = CO3-2 + Ca+2"                    },
        { "Strontianite", "Strontianite = Sr+2 + CO3-2"                 },
        { "Witherite"   , "Witherite = Ba+2 + CO3-2"                    },
        { "Gypsum"      , "Gypsum = Ca+2 + SO4-2 + 2*H2O"               },
        { "Anhydrite"   , "Anhydrite = Ca+2 + SO4-2"                    },
        { "Celestite"   , "Celestite = Sr+2 + SO4-2"                    },
        { "Barite"      , "Barite = Ba+2 + SO4-2"                       },
        { "Fluorite"    , "Fluorite = Ca+2 + 2*F-"                      },
        { "SiO2(a)"     , "SiO2(a) + 2*H2O = H4SiO4"                    },
        { "Chalcedony"  , "Chalcedony + 2*H2O = H4SiO4"                 },
        { "Quartz"      , "Quartz + 2*H2O = H4SiO4"                     },
        { "Chrysotile"  , "Chrysotile + 6*H+ = H2O + 2*H4SiO4 + 3*Mg+2" },
        { "CO2(g)"      , "CO2(g) = CO2"                                },
        { "H2O(g)"      , "H2O(g) = H2O"                                },
        { "O2(g)"       , "O2(g) = O2"                                  },
        { "H2(g)"       , "H2(g) = H2"                                  },
        { "N2(g)"       , "N2(g) = N2"                                  },
        { "H2S(g)"      , "H2S(g) =  H+ + HS-"                          },
        { "CH4(g)"      , "CH4(g) = CH4"                                },
        { "NH3(g)"      , "NH3(g) = NH3"                                },
        { "Oxg(g)"      , "Oxg(g) = Oxg"                                },
        { "Hdg(g)"      , "Hdg(g) = Hdg"                                },
        { "Ntg(g)"      , "Ntg(g) = Ntg"                                },
        { "Mtg(g)"      , "Mtg(g) = Mtg"                                },
        { "H2Sg(g)"     , "H2Sg(g) =  H+ + HSg-"                        },
        { "Melanterite" , "Melanterite = 7*H2O + Fe+2 + SO4-2"          },
    };

    for(auto [species, reactionstr] : reactions)
    {
        const auto reaction = parseReactionEquation(reactionstr);
        INFO("species: " << species)
        INFO("T: " << T << " K")
        INFO("P: " << P << " Pa")
        CHECK( lgK(db, T, P, species) == Approx(lgK_fromG0(db, T, P, reaction) * lgK_sign(db, species)).scale(1.0) );
    }
}

TEST_CASE("Testing pressure correction in standard thermodynamic properties calculations", "[PhreeqcDatabase]")
{
    static const auto db = test::getPhreeqcDatabase("phreeqc.dat");

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

        INFO("species: " << species.name())

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

        INFO("species: " << species.name())

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

        INFO("species: " << species.name())

        CHECK(propsT0P0.G0  == Approx(-95213.7) );
        CHECK(propsT0P0.H0  == Approx(-24010.3) );
        CHECK(propsT0P0.V0  == Approx(3.44329e-05) );
        CHECK(propsT0P0.Cp0 == Approx(657.168) );

        CHECK(propsT1P0.G0  == Approx(-113039) );
        CHECK(propsT1P0.H0  == Approx(12124.6) );
        CHECK(propsT1P0.V0  == Approx(3.77425e-05) );
        CHECK(propsT1P0.Cp0 == Approx(574.514) );

        CHECK(propsT0P1.G0  == Approx(-95161) );
        CHECK(propsT0P1.H0  == Approx(-23957.6) );
        CHECK(propsT0P1.V0  == Approx(3.44286e-05) );
        CHECK(propsT0P1.Cp0 == Approx(657.168) );

        CHECK(propsT1P1.G0  == Approx(-112978) );
        CHECK(propsT1P1.H0  == Approx(12185.4) );
        CHECK(propsT1P1.V0  == Approx(3.77222e-05) );
        CHECK(propsT1P1.Cp0 == Approx(574.514) );
    }

    {
        const auto species = db.species().get("Calcite");

        const auto propsT0P0 = species.props(T0, P0);
        const auto propsT1P0 = species.props(T1, P0);
        const auto propsT0P1 = species.props(T0, P1);
        const auto propsT1P1 = species.props(T1, P1);

        INFO("species: " << species.name())

        CHECK(propsT0P0.G0  == Approx(-48403.7) );
        CHECK(propsT0P0.H0  == Approx(12253.3) );
        CHECK(propsT0P0.V0  == Approx(3.69e-05) );
        CHECK(propsT0P0.Cp0 == Approx(531.164) );

        CHECK(propsT1P0.G0  == Approx(-63817.2) );
        CHECK(propsT1P0.H0  == Approx(47329.9) );
        CHECK(propsT1P0.V0  == Approx(3.69e-05) );
        CHECK(propsT1P0.Cp0 == Approx(638.057) );

        CHECK(propsT0P1.G0  == Approx(-48348.6) );
        CHECK(propsT0P1.H0  == Approx(12308.4) );
        CHECK(propsT0P1.V0  == Approx(3.69e-05) );
        CHECK(propsT0P1.Cp0 == Approx(531.164) );

        CHECK(propsT1P1.G0  == Approx(-63756.5) );
        CHECK(propsT1P1.H0  == Approx(47390.5) );
        CHECK(propsT1P1.V0  == Approx(3.69e-05) );
        CHECK(propsT1P1.Cp0 == Approx(638.057) );
    }
}

//=================================================================================================
// AUXILIARY FUNCTIONS: DEFINITION
//=================================================================================================

auto evalReactionStandardThermoProps(const Species& species, real T, real P) -> ReactionStandardThermoProps
{
    const auto dV0 = 0.0;
    return species.reaction().reactionThermoModel()({T, P, dV0});
}

auto lgK(const PhreeqcDatabase& db, real T, real P, const String& name) -> real
{
    const auto species = db.species().get(name);
    const auto dG0 = evalReactionStandardThermoProps(species, T, P).dG0;
    const auto R = universalGasConstant;
    const auto lnK = -dG0/(R*T);
    return lnK / ln10;
}

auto dH0(const PhreeqcDatabase& db, real T, real P, const String& name) -> real
{
    const auto species = db.species().get(name);
    const auto dH0 = evalReactionStandardThermoProps(species, T, P).dH0;
    return dH0;
}

auto dCp0(const PhreeqcDatabase& db, real T, real P, const String& name) -> real
{
    const auto species = db.species().get(name);
    const auto dCp0 = evalReactionStandardThermoProps(species, T, P).dCp0;
    return dCp0;
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

auto lgK_analytic(real T, Vec<real> A) -> real
{
    assert(A.size() <= 6);
    for(auto i = A.size(); i < 6; ++i) // A must have 6 entries -- fill out missing ones with zero!
        A.push_back(0.0);
    return A[0] + A[1]*T + A[2]/T + A[3]*log10(T) + A[4]/(T*T) + A[5]*T*T;
}

auto lgK_sign(const PhreeqcDatabase& db, const String& name) -> real
{
    // NOTE: The logic below is needed because in Reaktoro any reaction
    // associated with a chemical species is represented as a formation
    // reaction, i.e., a reaction in which the species itself is on the
    // product side (right-hand side of the reaction). In PHREEQC, this does
    // not for gases and minerals. That's why we multiply the log(K) by -1 for
    // chemical species with gas and solid aggregate state.
    const auto species = db.species().get(name);
    const auto is_phreeqc_phase = oneof(species.aggregateState(), AggregateState::Gas, AggregateState::Solid);
    const auto sign = is_phreeqc_phase ? -1.0 : 1.0;
    return sign;
}

auto lgK_fromG0(const PhreeqcDatabase& db, real T, real P, const Pairs<String, double>& reaction) -> real
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
    return R"#(
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
    )#";
}
