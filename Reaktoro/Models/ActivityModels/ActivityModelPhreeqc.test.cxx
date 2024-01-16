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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Warnings.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPhreeqc.hpp>
using namespace Reaktoro;

// -----------------------------------------------------------------------------------------------------------------
// NOTE: For the three tests in this file, we used the following PHREEQC
// scripts, with the databases: phreeqc.dat, lnll.dat, and thermoddem-v1.10.dat,
// respectively. The scripts were executed using `phreeqc script.in script.out
// databasefile`. The data was collected from the selected.out file using Visual
// Studio Code. Zero molalities reported by PHREEQC were replaced by 1e-99. This
// is because the activity model implemented in Reaktoro does not accept zero
// amounts and zero mole fractions. The molalities were assumed mole amounts in
// the tests below because PHREEQC maintained at all cases 1 kg of water in the
// solution.
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SOLUTION
//     temperature 60.0
//     pressure 100.0
//     units mol/kgw
//     pH 7.0 charge
//     B  0.10
//     Ba 1.20
//     Br 0.20
//     C  4.00
//     Ca 0.50
//     Cl 2.00
//     Fe 0.70
//     K  0.80
//     Li 0.60
//     Mg 0.40
//     Mn 0.60
//     Na 2.00
//     S  0.50
//     Si 0.10
//     Sr 1.40
// SELECTED_OUTPUT
//     -file selected.out
//     -high_precision true
//     -molalities OH- H+ H2O BH4- BaB(OH)4+ BO2- NaB(OH)4 CaB(OH)4+ MgB(OH)4+ B(OH)3 B2O(OH)5- Ba+2 BaCO3 BaOH+ BaCl+ Br3- Br- NaBr KBr Br2 BrO- HBrO BrO3- BrO4- C2H4 C2H6 CH4 CO SrCO3 CaCO3 CO3-2 NaCO3- MgCO3 NaHCO3 HCO3- CaHCO3+ MgHCO3+ CO2 MnCO3 FeCO3 MnHCO3+ FeCO3+ FeHCO3+ Ca+2 CaOH+ CaSO4 CaCl2 CaCl+ Cl- NaCl SrCl+ KCl LiCl MgCl+ HCl MnCl+ MnCl3- FeCl+ FeCl4-2 FeCl2 FeCl2+ FeCl+2 FeCl4- ClO- HClO ClO2- HClO2 ClO3- ClO4- Fe(OH)3- Fe(OH)4-2 Fe(OH)2 FeOH+ Fe+2 FeSO4 Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 Fe+3 FeSO4+ Fe(SO4)2- Fe2(OH)2+4 Fe3(OH)4+5 H2 K+ KOH KSO4- KHSO4 Li+ LiOH LiSO4- Mg4(OH)4+4 MgSO4 Mg+2 Mn(OH)4-2 Mn(OH)3- Mn(OH)2 MnOH+ MnSO4 Mn+2 Mn2(OH)3+ Mn2OH+3 Mn+3 MnO4-2 MnO4- Na+ NaOH NaHSiO3 NaSO4- O2 S-2 HS- H2S S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- S2O4-2 SO3-2 HSO3- SO2 H2SO3 S2O6-2 S3O6-2 S4O6-2 S5O6-2 S2O5-2 SrSO4 SO4-2 HSO4- H2SO4 S2O8-2 HSO5- H2SiO4-2 HSiO3- H4(H2SiO4)4-4 SiO2 H6(H2SiO4)4-2 Sr+2 SrOH+
//     -activities OH- H+ H2O BH4- BaB(OH)4+ BO2- NaB(OH)4 CaB(OH)4+ MgB(OH)4+ B(OH)3 B2O(OH)5- Ba+2 BaCO3 BaOH+ BaCl+ Br3- Br- NaBr KBr Br2 BrO- HBrO BrO3- BrO4- C2H4 C2H6 CH4 CO SrCO3 CaCO3 CO3-2 NaCO3- MgCO3 NaHCO3 HCO3- CaHCO3+ MgHCO3+ CO2 MnCO3 FeCO3 MnHCO3+ FeCO3+ FeHCO3+ Ca+2 CaOH+ CaSO4 CaCl2 CaCl+ Cl- NaCl SrCl+ KCl LiCl MgCl+ HCl MnCl+ MnCl3- FeCl+ FeCl4-2 FeCl2 FeCl2+ FeCl+2 FeCl4- ClO- HClO ClO2- HClO2 ClO3- ClO4- Fe(OH)3- Fe(OH)4-2 Fe(OH)2 FeOH+ Fe+2 FeSO4 Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 Fe+3 FeSO4+ Fe(SO4)2- Fe2(OH)2+4 Fe3(OH)4+5 H2 K+ KOH KSO4- KHSO4 Li+ LiOH LiSO4- Mg4(OH)4+4 MgSO4 Mg+2 Mn(OH)4-2 Mn(OH)3- Mn(OH)2 MnOH+ MnSO4 Mn+2 Mn2(OH)3+ Mn2OH+3 Mn+3 MnO4-2 MnO4- Na+ NaOH NaHSiO3 NaSO4- O2 S-2 HS- H2S S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- S2O4-2 SO3-2 HSO3- SO2 H2SO3 S2O6-2 S3O6-2 S4O6-2 S5O6-2 S2O5-2 SrSO4 SO4-2 HSO4- H2SO4 S2O8-2 HSO5- H2SiO4-2 HSiO3- H4(H2SiO4)4-4 SiO2 H6(H2SiO4)4-2 Sr+2 SrOH+
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SOLUTION
//     temperature 60.0
//     pressure 100.0
//     units mol/kgw
//     pH 7.0 charge
//     B  0.10
//     Ba 1.20
//     Br 0.20
//     C  4.00
//     Ca 0.50
//     Cl 2.00
//     Fe 0.70
//     K  0.80
//     Li 0.60
//     Mg 0.40
//     Mn 0.60
//     Na 2.00
//     S  0.50
//     Si 0.10
//     Sr 1.40
// SELECTED_OUTPUT
//     -file selected.out
//     -high_precision true
//     -molalities OH- H+ H2O BH4- BaB(OH)4+ BO2- NaB(OH)4 CaB(OH)4+ MgB(OH)4+ B(OH)3 B2O(OH)5- Ba+2 BaCO3 BaOH+ BaCl+ Br3- Br- NaBr KBr Br2 BrO- HBrO BrO3- BrO4- C2H4 C2H6 CH4 CO SrCO3 CaCO3 CO3-2 NaCO3- MgCO3 NaHCO3 HCO3- CaHCO3+ MgHCO3+ CO2 MnCO3 FeCO3 MnHCO3+ FeCO3+ FeHCO3+ Ca+2 CaOH+ CaSO4 CaCl2 CaCl+ Cl- NaCl SrCl+ KCl LiCl MgCl+ HCl MnCl+ MnCl3- FeCl+ FeCl4-2 FeCl2 FeCl2+ FeCl+2 FeCl4- ClO- HClO ClO2- HClO2 ClO3- ClO4- Fe(OH)3- Fe(OH)4-2 Fe(OH)2 FeOH+ Fe+2 FeSO4 Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 Fe+3 FeSO4+ Fe(SO4)2- Fe2(OH)2+4 Fe3(OH)4+5 H2 K+ KOH KSO4- KHSO4 Li+ LiOH LiSO4- Mg4(OH)4+4 MgSO4 Mg+2 Mn(OH)4-2 Mn(OH)3- Mn(OH)2 MnOH+ MnSO4 Mn+2 Mn2(OH)3+ Mn2OH+3 Mn+3 MnO4-2 MnO4- Na+ NaOH NaHSiO3 NaSO4- O2 S-2 HS- H2S S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- S2O4-2 SO3-2 HSO3- SO2 H2SO3 S2O6-2 S3O6-2 S4O6-2 S5O6-2 S2O5-2 SrSO4 SO4-2 HSO4- H2SO4 S2O8-2 HSO5- H2SiO4-2 HSiO3- H4(H2SiO4)4-4 SiO2 H6(H2SiO4)4-2 Sr+2 SrOH+
//     -activities OH- H+ H2O BH4- BaB(OH)4+ BO2- NaB(OH)4 CaB(OH)4+ MgB(OH)4+ B(OH)3 B2O(OH)5- Ba+2 BaCO3 BaOH+ BaCl+ Br3- Br- NaBr KBr Br2 BrO- HBrO BrO3- BrO4- C2H4 C2H6 CH4 CO SrCO3 CaCO3 CO3-2 NaCO3- MgCO3 NaHCO3 HCO3- CaHCO3+ MgHCO3+ CO2 MnCO3 FeCO3 MnHCO3+ FeCO3+ FeHCO3+ Ca+2 CaOH+ CaSO4 CaCl2 CaCl+ Cl- NaCl SrCl+ KCl LiCl MgCl+ HCl MnCl+ MnCl3- FeCl+ FeCl4-2 FeCl2 FeCl2+ FeCl+2 FeCl4- ClO- HClO ClO2- HClO2 ClO3- ClO4- Fe(OH)3- Fe(OH)4-2 Fe(OH)2 FeOH+ Fe+2 FeSO4 Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 Fe+3 FeSO4+ Fe(SO4)2- Fe2(OH)2+4 Fe3(OH)4+5 H2 K+ KOH KSO4- KHSO4 Li+ LiOH LiSO4- Mg4(OH)4+4 MgSO4 Mg+2 Mn(OH)4-2 Mn(OH)3- Mn(OH)2 MnOH+ MnSO4 Mn+2 Mn2(OH)3+ Mn2OH+3 Mn+3 MnO4-2 MnO4- Na+ NaOH NaHSiO3 NaSO4- O2 S-2 HS- H2S S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- S2O4-2 SO3-2 HSO3- SO2 H2SO3 S2O6-2 S3O6-2 S4O6-2 S5O6-2 S2O5-2 SrSO4 SO4-2 HSO4- H2SO4 S2O8-2 HSO5- H2SiO4-2 HSiO3- H4(H2SiO4)4-4 SiO2 H6(H2SiO4)4-2 Sr+2 SrOH+
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SOLUTION
//     temperature 60.0
//     pressure 100.0
//     units mol/kgw
//     pH 7.0 charge
//     B  0.10
//     Ba 1.20
//     Br 0.20
//     C  4.00
//     Ca 0.50
//     Cl 2.00
//     Fe 0.70
//     K  0.80
//     Li 0.60
//     Mg 0.40
//     Mn 0.60
//     Na 2.00
//     S  0.50
//     Si 0.10
//     Sr 1.40
// SELECTED_OUTPUT
//     -file selected.out
//     -high_precision true
//     -molalities OH- H+ H2O B(OH)3 B(OH)4- NaB(OH)4 Ba+2 BaCO3 Ba(HCO3)+ BaCl+ BaOH+ Br3- Br- NaBr KBr HBr BrO- HBrO BrO3- BrO4- CH4 CO CaCO3 MgCO3 HCO3- NaCO3- NaHCO3 Mg(HCO3)+ CO3-2 Ca(HCO3)+ Sr(CO3) Sr(HCO3)+ BaCO3 Ba(HCO3)+ Fe(CO3)2- CO2 FeCO3 Fe(CO3)2-2 FeCO3OH- FeHCO3+ FeCO3+ CaCO3 Ca+2 CaSO4 Ca(HCO3)+ CaCl+ CaCl2 CaOH+ Cl- KCl MgCl+ MnCl+ CaCl+ BaCl+ SrCl+ LiCl CaCl2 FeCl+ HCl FeCl2 FeCl+2 FeCl2+ ClO- HClO ClO2- HClO2 ClO2 ClO3- ClO4- FeCO3 Fe(CO3)2-2 FeCO3OH- FeSO4 Fe+2 FeHCO3+ FeOH+ FeCl+ FeO FeCl2 HFeO2- FeHSO4+ Fe(HS)2 Fe(OH)4- Fe(CO3)2- HFeO2 FeO+ FeOH+2 FeCO3+ FeSO4+ FeCl+2 Fe+3 Fe2(OH)2+4 FeCl2+ FeHSO4+2 H2 K+ KCl KSO4- KBr KOH Li+ LiCl LiOH MgCO3 Mg+2 MgSO4 Mg(HCO3)+ MgCl+ MgOH+ Mg4(OH)4+4 Mn+2 MnSO4 MnCl+ MnOH+ MnO HMnO2- MnO2-2 Mn+3 MnO4-2 MnO4- Na+ NaCO3- NaSO4- NaHCO3 NaB(OH)4 NaBr NaOH O2 HS- H2S S-2 S2-2 Fe(HS)2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- H2S2O3 S2-2 S4O6-2 S3-2 S5O6-2 S4-2 S5-2 S2O4-2 HS2O4- H2S2O4 S3O6-2 S4O6-2 SO3-2 HSO3- H2SO3 S2O5-2 S3O6-2 S2O6-2 SO4-2 NaSO4- MgSO4 CaSO4 KSO4- MnSO4 SrSO4 FeSO4 HSO4- FeSO4+ FeHSO4+ FeHSO4+2 S2O8-2 HSO5- H4SiO4 HSiO3- Si4O6(OH)6-2 Si2O2(OH)5- Si2O3(OH)4-2 Si3O6(OH)3-3 Si3O5(OH)5-3 H2SiO4-2 Si6O15-6 Si4O8(OH)4-4 Si4O7(OH)6-4 Si4O12H4-4 Sr(CO3) Sr(HCO3)+ Sr+2 SrSO4 SrCl+ SrOH+
//     -activities OH- H+ H2O B(OH)3 B(OH)4- NaB(OH)4 Ba+2 BaCO3 Ba(HCO3)+ BaCl+ BaOH+ Br3- Br- NaBr KBr HBr BrO- HBrO BrO3- BrO4- CH4 CO CaCO3 MgCO3 HCO3- NaCO3- NaHCO3 Mg(HCO3)+ CO3-2 Ca(HCO3)+ Sr(CO3) Sr(HCO3)+ BaCO3 Ba(HCO3)+ Fe(CO3)2- CO2 FeCO3 Fe(CO3)2-2 FeCO3OH- FeHCO3+ FeCO3+ CaCO3 Ca+2 CaSO4 Ca(HCO3)+ CaCl+ CaCl2 CaOH+ Cl- KCl MgCl+ MnCl+ CaCl+ BaCl+ SrCl+ LiCl CaCl2 FeCl+ HCl FeCl2 FeCl+2 FeCl2+ ClO- HClO ClO2- HClO2 ClO2 ClO3- ClO4- FeCO3 Fe(CO3)2-2 FeCO3OH- FeSO4 Fe+2 FeHCO3+ FeOH+ FeCl+ FeO FeCl2 HFeO2- FeHSO4+ Fe(HS)2 Fe(OH)4- Fe(CO3)2- HFeO2 FeO+ FeOH+2 FeCO3+ FeSO4+ FeCl+2 Fe+3 Fe2(OH)2+4 FeCl2+ FeHSO4+2 H2 K+ KCl KSO4- KBr KOH Li+ LiCl LiOH MgCO3 Mg+2 MgSO4 Mg(HCO3)+ MgCl+ MgOH+ Mg4(OH)4+4 Mn+2 MnSO4 MnCl+ MnOH+ MnO HMnO2- MnO2-2 Mn+3 MnO4-2 MnO4- Na+ NaCO3- NaSO4- NaHCO3 NaB(OH)4 NaBr NaOH O2 HS- H2S S-2 S2-2 Fe(HS)2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- H2S2O3 S2-2 S4O6-2 S3-2 S5O6-2 S4-2 S5-2 S2O4-2 HS2O4- H2S2O4 S3O6-2 S4O6-2 SO3-2 HSO3- H2SO3 S2O5-2 S3O6-2 S2O6-2 SO4-2 NaSO4- MgSO4 CaSO4 KSO4- MnSO4 SrSO4 FeSO4 HSO4- FeSO4+ FeHSO4+ FeHSO4+2 S2O8-2 HSO5- H4SiO4 HSiO3- Si4O6(OH)6-2 Si2O2(OH)5- Si2O3(OH)4-2 Si3O6(OH)3-3 Si3O5(OH)5-3 H2SiO4-2 Si6O15-6 Si4O8(OH)4-4 Si4O7(OH)6-4 Si4O12H4-4 Sr(CO3) Sr(HCO3)+ Sr+2 SrSO4 SrCl+ SrOH+
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/// The data used to test ActivityModelPhreeqc.
struct TestDataActivityModelPhreeqc
{
    /// The name of the PHREEQC database that should be used.
    String dbname;

    /// The list of aqueous species in the aqueous solution.
    StringList species;

    /// The temperature used to produce the benchmark data through PHREEQC script (in °C)
    real TdegC;

    /// The pressure used to produce the benchmark data through PHREEQC script (in atn)
    real Patm;

    /// The amounts of the species in the solution (in mol)
    ArrayXr n;

    /// The activities of the species in the solution (in log base 10) obtained from PHREEQC.
    ArrayXr la;
};

TEST_CASE("Testing ActivityModelPhreeqc", "[ActivityModelPhreeqc]")
{
    TestDataActivityModelPhreeqc data1;
    data1.dbname  = "phreeqc.dat";
    data1.species = "OH- H+ H2O H3BO3 H2BO3- Ba+2 BaSO4 BaCO3 BaHCO3+ BaOH+ Br- CH4 CaCO3 MnCO3 NaCO3- HCO3- MgCO3 MgHCO3+ CaHCO3+ SrHCO3+ CO3-2 SrCO3 MnHCO3+ NaHCO3 CO2 FeCO3 FeHCO3+ (CO2)2 Ca+2 CaSO4 CaOH+ CaHSO4+ Cl- MnCl+ MnCl2 MnCl3- FeCl+ FeCl+2 FeCl2+ FeCl3 Fe+2 FeOH+ FeSO4 Fe(OH)2 Fe(OH)3- FeHSO4+ Fe(HS)2 Fe(HS)3- Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 FeSO4+ Fe(SO4)2- Fe+3 Fe2(OH)2+4 Fe3(OH)4+5 FeHSO4+2 H2 K+ KSO4- Li+ LiSO4- MgSO4 Mg+2 MgOH+ Mn+2 MnSO4 MnOH+ Mn(OH)3- Mn+3 Na+ NaSO4- NaOH O2 HS- S-2 H2S (H2S)2 SO4-2 SrSO4 HSO4- H4SiO4 H3SiO4- H2SiO4-2 Sr+2 SrOH+";
    data1.TdegC   = 60.0;
    data1.Patm    = 100.0;
    data1.n       = ArrayXr{{ 1.045246946269e-04, 2.569031392499e-09, 5.550621669627e+01, 5.117559439648e-02, 4.882440560352e-02, 4.566925050950e-02, 2.813774490079e-02, 1.315144774072e-02, 1.304127753295e-02, 2.793160509090e-07, 1.000000000000e-01, 1.000000000000e-99, 2.499216738509e-01, 1.638690339138e-01, 1.203657410487e-01, 9.700728350679e-02, 9.628983342313e-02, 7.441337689308e-02, 6.617950814229e-02, 3.999103417931e-02, 2.044288187112e-02, 1.867830652383e-02, 1.508182043098e-02, 1.145575741425e-02, 1.107215719740e-04, 2.595563477225e-07, 4.022883325544e-08, 1.085471911400e-09, 1.303740721553e-01, 5.352203284241e-02, 2.711205152830e-06, 1.803871250741e-09, 9.922306790609e-01, 6.275508147571e-03, 4.131358217017e-04, 2.225120176913e-04, 5.095069439034e-09, 1.829951796501e-14, 3.482590616919e-15, 6.412737755430e-17, 5.693605526519e-08, 2.408882682313e-08, 1.606321925517e-08, 5.184280305489e-10, 4.561040614460e-11, 3.670562378786e-16, 1.000000000000e-99, 1.000000000000e-99, 4.214022180612e-01, 7.769670962179e-02, 9.006615316525e-04, 8.252894766161e-09, 1.886194001125e-14, 2.506300742958e-15, 9.412197371089e-16, 1.286167922323e-19, 9.029580342171e-23, 2.366932808311e-23, 8.253943747732e-30, 5.356678064722e-01, 6.433219352777e-02, 9.474572460148e-02, 5.254275398515e-03, 1.313525238432e-01, 9.584938762260e-02, 2.094878217999e-03, 1.073407670050e-02, 2.942077032438e-03, 4.618359289058e-04, 6.398184726758e-12, 1.012250466290e-21, 7.908700133118e-01, 7.730848822529e-02, 1.953161083459e-15, 8.281023564144e-26, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.253765391762e-01, 1.177410386960e-02, 3.316687643813e-09, 5.464864351916e-02, 4.530089932422e-02, 5.045715661899e-05, 2.955621076764e-02, 3.446596178519e-07 }};
    data1.la      = ArrayXr{{ -4.286399139975e+00, -8.740396300458e+00, -3.886360917357e-02, -1.018722512235e+00, -1.206042449304e+00, -2.041209329790e+00, -1.278496119893e+00, -1.608811844297e+00, -1.779359276842e+00, -6.792234610668e+00, -1.337367563669e+00, -8.565298326199e+01, -3.299814869654e-01, -5.132885145155e-01, -8.141765201461e-01, -1.238309359906e+00, -7.442049724631e-01, -1.407679371676e+00, -1.387102444607e+00, -1.623151067866e+00, -2.589912693407e+00, -1.456447909665e+00, -2.059876738205e+00, -1.668761599368e+00, -3.683553164791e+00, -6.313553753389e+00, -7.290141977078e+00, -8.692166818100e+00, -1.388976707656e+00, -9.992528076697e-01, -5.461517032886e+00, -8.638473877085e+00, -2.554295096079e-01, -2.440681606075e+00, -3.111692554949e+00, -3.890977031999e+00, -8.187529306882e+00, -1.469088236540e+01, -1.469642807826e+01, -1.592074192776e+01, -8.075916544391e+00, -7.856514854385e+00, -7.521952820551e+00, -9.013096932777e+00, -1.057926656434e+01, -1.532994680510e+01, -1.662132988390e+02, -2.477063405298e+02, -6.004168859236e-01, -8.373827806087e-01, -3.270552089226e+00, -9.036715707956e+00, -1.396274414455e+01, -1.449564623127e+01, -1.637780037927e+01, -1.720557295209e+01, -1.941131778649e+01, -2.220453172699e+01, -2.881112390311e+01, -5.358903124236e-01, -1.416685343377e+00, -1.231266358136e+00, -2.517817671603e+00, -6.093469860058e-01, -1.362123141480e+00, -2.874167073445e+00, -2.800539224540e+00, -2.259131368183e+00, -3.573842787646e+00, -1.143227372907e+01, -2.234620341941e+01, -1.542172523615e-01, -1.336886522437e+00, -1.443704734544e+01, -2.480970138741e+01, -8.356648158620e+01, -8.677587654272e+01, -8.571339100274e+01, -1.723369020101e+02, -1.963991574559e+00, -1.656857545079e+00, -8.373974841293e+00, -9.902060214368e-01, -1.623223561210e+00, -5.197532039145e+00, -2.118249737552e+00, -6.700940101693e+00 }};

    TestDataActivityModelPhreeqc data2;
    data2.dbname  = "llnl.dat";
    data2.species = "OH- H+ H2O BH4- BaB(OH)4+ BO2- NaB(OH)4 CaB(OH)4+ MgB(OH)4+ B(OH)3 B2O(OH)5- Ba+2 BaCO3 BaOH+ BaCl+ Br3- Br- NaBr KBr Br2 BrO- HBrO BrO3- BrO4- C2H4 C2H6 CH4 CO SrCO3 CaCO3 CO3-2 NaCO3- MgCO3 NaHCO3 HCO3- CaHCO3+ MgHCO3+ CO2 MnCO3 FeCO3 MnHCO3+ FeCO3+ FeHCO3+ Ca+2 CaOH+ CaSO4 CaCl2 CaCl+ Cl- NaCl SrCl+ KCl LiCl MgCl+ HCl MnCl+ MnCl3- FeCl+ FeCl4-2 FeCl2 FeCl2+ FeCl+2 FeCl4- ClO- HClO ClO2- HClO2 ClO3- ClO4- Fe(OH)3- Fe(OH)4-2 Fe(OH)2 FeOH+ Fe+2 FeSO4 Fe(OH)4- Fe(OH)3 Fe(OH)2+ FeOH+2 Fe+3 FeSO4+ Fe(SO4)2- Fe2(OH)2+4 Fe3(OH)4+5 H2 K+ KOH KSO4- KHSO4 Li+ LiOH LiSO4- Mg4(OH)4+4 MgSO4 Mg+2 Mn(OH)4-2 Mn(OH)3- Mn(OH)2 MnOH+ MnSO4 Mn+2 Mn2(OH)3+ Mn2OH+3 Mn+3 MnO4-2 MnO4- Na+ NaOH NaHSiO3 NaSO4- O2 S-2 HS- H2S S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- S2O4-2 SO3-2 HSO3- SO2 H2SO3 S2O6-2 S3O6-2 S4O6-2 S5O6-2 S2O5-2 SrSO4 SO4-2 HSO4- H2SO4 S2O8-2 HSO5- H2SiO4-2 HSiO3- H4(H2SiO4)4-4 SiO2 H6(H2SiO4)4-2 Sr+2 SrOH+";
    data2.TdegC   = 60.0;
    data2.Patm    = 100.0;
    data2.n       = ArrayXr{{ 1.918731734567e-04, 7.782597999947e-10, 5.552532510078e+01, 1.000000000000e-99, 3.281285670687e-03, 3.671808336828e-02, 1.105382262851e-02, 1.362984270422e-02, 1.976515385560e-02, 1.555181177249e-02, 1.073007339825e-13, 4.641865420494e-02, 4.708253499380e-02, 4.523379763063e-07, 3.217072792597e-03, 1.753084701488e-28, 9.733189377857e-02, 2.025481255730e-03, 6.426249656976e-04, 4.409290465748e-28, 2.019166053607e-24, 2.301077924261e-25, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 3.035821582364e-35, 5.373135712576e-02, 3.709406914142e-01, 5.139034798827e-02, 5.552854819577e-03, 1.532088849042e-01, 2.373910846380e-02, 9.098909506892e-02, 1.580411410079e-02, 1.943496105986e-02, 5.653365876692e-05, 1.672043008264e-01, 1.330667518702e-05, 8.490363944095e-04, 7.441436228186e-10, 2.871761829621e-06, 6.152672450377e-02, 3.272727551718e-06, 3.482530050212e-02, 8.115572259980e-04, 2.458496821342e-03, 9.000053991978e-01, 6.180839884659e-02, 3.193243738651e-03, 1.027130373107e-02, 1.604544532324e-03, 1.025198104296e-02, 8.105131054048e-11, 5.182462048039e-03, 1.279774782494e-04, 4.481812712322e-08, 1.311471210644e-09, 1.083945172262e-10, 8.238894525012e-15, 5.641469101035e-16, 2.802378614452e-18, 4.659469586880e-27, 8.042587066048e-29, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 2.069850295760e-11, 1.800534891845e-16, 2.637580907743e-10, 3.748799542230e-08, 3.148082903232e-07, 1.422653628356e-07, 1.653120497998e-01, 3.338653632946e-01, 8.058566487950e-04, 9.989900132505e-09, 2.507998968477e-15, 1.310130047425e-15, 7.866676329233e-17, 1.010232105839e-13, 4.774566051834e-13, 1.217470952978e-29, 5.203002801876e-01, 1.408838584300e-06, 6.878438217824e-02, 9.877700253186e-11, 8.921635121997e-02, 2.118719756070e-06, 9.176985527954e-03, 5.620947856477e-08, 1.459155102725e-01, 5.142328402693e-02, 1.838990263979e-14, 2.496415232905e-10, 1.350155154190e-07, 6.209692667032e-05, 7.131345379763e-03, 6.415404566821e-03, 6.508093520628e-03, 5.527036592613e-06, 3.503073887546e-22, 1.721126753270e-20, 7.124769586109e-25, 7.558192138820e-01, 1.078659192784e-05, 6.545571535101e-02, 7.453461816088e-02, 1.792688255154e-25, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 4.909267289015e-31, 1.873449739676e-33, 1.729581764823e-40, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.304684450915e-02, 1.465848667868e-01, 4.318440097529e-09, 6.857696416003e-22, 1.000000000000e-99, 1.926123022732e-38, 2.064574111479e-06, 7.275471799705e-03, 6.651323969222e-04, 1.098081882101e-02, 3.406349966619e-03, 3.002811173374e-02, 4.428927069978e-07 }};
    data2.la      = ArrayXr{{ -3.915556115534e+00, -9.158430923754e+00, -3.740756498489e-02, -1.551317902964e+02, -2.657665744536e+00, -1.608829782938e+00, -1.956487508460e+00, -2.039218942635e+00, -1.877809587184e+00, -1.808219008762e+00, -1.314310709377e+01, -2.169270094543e+00, -1.327140162137e+00, -6.518246735774e+00, -2.666248898725e+00, -2.792990688667e+01, -1.240161007768e+00, -2.693471771523e+00, -3.192042406307e+00, -2.735563129071e+01, -2.386853775028e+01, -2.463806867403e+01, -5.397304915097e+01, -8.064023126233e+01, -1.732245711341e+02, -1.592745841216e+02, -8.932006369925e+01, -3.451772375586e+01, -1.269772190051e+00, -4.306955228663e-01, -2.197154939988e+00, -2.429193467594e+00, -8.147160483262e-01, -1.624535595260e+00, -1.214720440724e+00, -1.974939630013e+00, -1.885126111694e+00, -4.059707665305e+00, -7.767525558374e-01, -4.875930444023e+00, -3.244783479591e+00, -9.302053022225e+00, -5.715561367776e+00, -1.929823292642e+00, -5.658799933873e+00, -1.458105127428e+00, -3.090680850996e+00, -2.783040135389e+00, -2.741710664460e-01, -1.208952506709e+00, -2.669477717226e+00, -1.988374428139e+00, -2.794648225048e+00, -2.162901992050e+00, -1.009123995853e+01, -2.459173656132e+00, -4.066576238205e+00, -7.522256082388e+00, -9.874390272367e+00, -9.964992684575e+00, -1.425784084348e+01, -1.615664428469e+01, -1.772618297638e+01, -2.650537330514e+01, -2.809460422889e+01, -4.593444313379e+01, -5.192307405755e+01, -5.295624567287e+01, -6.405936495060e+01, -1.085777085075e+01, -1.673674749198e+01, -9.578794209515e+00, -7.599817568284e+00, -7.220840927052e+00, -6.846900824163e+00, -9.554052755172e-01, -4.764286342859e-01, -3.267451993055e+00, -8.908475351823e+00, -1.583949871059e+01, -1.505639537934e+01, -1.627791850481e+01, -1.638695070365e+01, -1.733440269670e+01, -2.872655615011e+01, -5.121621214886e-01, -5.851138762720e+00, -1.336219945708e+00, -1.000534415699e+01, -1.154949842568e+00, -5.673926483800e+00, -2.211009739679e+00, -1.064156231369e+01, -8.358985417752e-01, -1.843913937191e+00, -1.472756960379e+01, -9.776392962561e+00, -6.869616321330e+00, -4.380639680098e+00, -2.146828529631e+00, -2.911663038867e+00, -2.360256001428e+00, -7.262302718966e+00, -2.346034576123e+01, -2.075633617844e+01, -2.434579955020e+01, -2.952918584672e-01, -4.967115750936e+00, -1.184052426469e+00, -1.301351755578e+00, -2.455850998531e+01, -8.684419210045e+01, -8.400267470690e+01, -8.650470177820e+01, -1.510194349984e+02, -2.151949272560e+02, -2.795682205217e+02, -3.441391245607e+02, -9.032135649226e+01, -9.846588741602e+01, -8.951413689778e+01, -3.121701982024e+01, -3.290106774004e+01, -3.976205890218e+01, -4.032068166775e+01, -5.433552201463e+01, -1.208671421592e+02, -1.725577259360e+02, -2.512820763384e+02, -7.069380585425e+01, -1.884494513458e+00, -1.826059897111e+00, -8.538382886767e+00, -2.116382174462e+01, -5.300312741277e+01, -3.788902576416e+01, -6.677318556539e+00, -2.311848624810e+00, -7.442997921171e+00, -1.959365274077e+00, -3.459859768678e+00, -2.358434616347e+00, -6.527411257578e+00 }};

    TestDataActivityModelPhreeqc data3;
    data3.dbname  = "thermoddem-v1.10.dat";
    data3.species = "OH- H+ H2O B(OH)3 B(OH)4- NaB(OH)4 Ba+2 BaCO3 Ba(HCO3)+ BaCl+ BaOH+ Br3- Br- NaBr KBr HBr BrO- HBrO BrO3- BrO4- CH4 CO CaCO3 MgCO3 HCO3- NaCO3- NaHCO3 Mg(HCO3)+ CO3-2 Ca(HCO3)+ Sr(CO3) Sr(HCO3)+ Fe(CO3)2- CO2 FeCO3 Fe(CO3)2-2 FeCO3OH- FeHCO3+ FeCO3+ Ca+2 CaSO4 CaCl+ CaCl2 CaOH+ Cl- KCl MgCl+ MnCl+ SrCl+ LiCl FeCl+ HCl FeCl2 FeCl+2 FeCl2+ ClO- HClO ClO2- HClO2 ClO2 ClO3- ClO4- FeSO4 Fe+2 FeOH+ FeO HFeO2- FeHSO4+ Fe(HS)2 Fe(OH)4- HFeO2 FeO+ FeOH+2 FeSO4+ Fe+3 Fe2(OH)2+4 FeHSO4+2 H2 K+ KSO4- KOH Li+ LiOH Mg+2 MgSO4 MgOH+ Mg4(OH)4+4 Mn+2 MnSO4 MnOH+ MnO HMnO2- MnO2-2 Mn+3 MnO4-2 MnO4- Na+ NaSO4- NaOH O2 HS- H2S S-2 S2-2 S3-2 S4-2 S5-2 S2O3-2 HS2O3- H2S2O3 S4O6-2 S5O6-2 S2O4-2 HS2O4- H2S2O4 S3O6-2 SO3-2 HSO3- H2SO3 S2O5-2 S2O6-2 SO4-2 SrSO4 HSO4- S2O8-2 HSO5- H4SiO4 HSiO3- Si4O6(OH)6-2 Si2O2(OH)5- Si2O3(OH)4-2 Si3O6(OH)3-3 Si3O5(OH)5-3 H2SiO4-2 Si6O15-6 Si4O8(OH)4-4 Si4O7(OH)6-4 Si4O12H4-4 Sr+2 SrOH+";
    data3.TdegC   = 60.0;
    data3.Patm    = 100.0;
    data3.n       = ArrayXr{{ 7.287078728321e-05, 1.867925061164e-09, 5.550929780738e+01, 4.538890839216e-02, 3.285361974499e-02, 2.175747186285e-02, 4.219441446692e-02, 2.807258376619e-02, 2.676211735246e-02, 2.965322034851e-03, 5.562379577486e-06, 1.025302933958e-27, 9.739113502515e-02, 2.000959695818e-03, 6.079052790305e-04, 6.894484306263e-18, 1.858323904866e-24, 4.112243542722e-25, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.468102361320e-34, 2.869280672732e-01, 1.439346731304e-01, 1.419837938187e-01, 1.254023883568e-01, 4.733295197580e-02, 4.413569626323e-02, 3.468496162469e-02, 2.999836882074e-02, 2.893599770138e-02, 2.840406684197e-02, 1.659684985161e-02, 2.200131860134e-04, 5.334303922442e-06, 1.943767863386e-06, 1.364801239554e-06, 3.344367620729e-08, 1.006669629514e-10, 9.724766129242e-02, 7.328164257675e-02, 1.128520900135e-02, 1.199564518830e-03, 5.948651668097e-05, 8.464523488320e-01, 6.251105543583e-02, 3.851848024526e-02, 3.152021348124e-02, 2.549360394147e-03, 1.798871000168e-03, 1.027371756252e-08, 1.292444344064e-10, 6.727941586468e-11, 1.083710967175e-14, 2.188870836113e-16, 4.218555997799e-27, 1.114368622869e-28, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 4.370373822319e-08, 3.550830947184e-08, 1.686113345811e-08, 7.918678543846e-10, 1.149815891012e-11, 2.751112406090e-17, 1.000000000000e-99, 4.685687832777e-01, 1.347841214249e-02, 1.347167933114e-03, 3.160132834067e-09, 4.150879521984e-14, 6.487494989689e-16, 3.172519095276e-16, 4.540269322875e-23, 1.396218689758e-29, 4.769021754053e-01, 5.997062129130e-02, 8.242588499129e-06, 9.819154420446e-02, 9.584795372373e-06, 9.811761447411e-02, 7.478671577181e-02, 4.923047197481e-04, 3.628848860290e-06, 1.118454293697e-01, 5.160748369896e-02, 4.956887614070e-03, 6.990315738873e-05, 8.255816795927e-08, 1.204723235527e-10, 7.307165189089e-21, 4.936882298898e-21, 5.720223768991e-25, 7.170112850207e-01, 8.649027155084e-02, 4.671537227467e-06, 1.160685314496e-25, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 1.000000000000e-99, 6.744261533533e-31, 4.251211164667e-33, 1.004358331586e-39, 1.000000000000e-99, 1.000000000000e-99, 1.390676965122e-01, 1.479551388693e-02, 1.100743995806e-08, 1.000000000000e-99, 4.601635902195e-38, 3.443163748302e-02, 1.407078542656e-02, 8.539217294132e-03, 7.235994103043e-03, 9.807382827629e-04, 1.464676759047e-04, 1.339393107078e-04, 8.552971399401e-06, 6.726080446733e-06, 1.568118954207e-06, 1.434007018025e-06, 1.276056118338e-06, 2.531062787178e-02, 4.433303799220e-06 }};
    data3.la      = ArrayXr{{ -4.322527873134e+00, -8.756197836487e+00, -3.780226955119e-02, -1.343050261958e+00, -1.662759415590e+00, -1.662391569443e+00, -2.212859377366e+00, -1.551717613594e+00, -1.726056516461e+00, -2.681505122449e+00, -5.408316365370e+00, -2.716849044167e+01, -1.228856919137e+00, -2.698761659020e+00, -3.216164085199e+00, -1.716149821264e+01, -2.391022122825e+01, -2.438592117268e+01, -5.417855936345e+01, -8.092267615075e+01, -8.840560271021e+01, -3.383324366281e+01, -5.422269670357e-01, -8.418345741293e-01, -9.967021000580e-01, -1.081036833857e+00, -1.324836409396e+00, -1.508787005025e+00, -2.374546093800e+00, -1.676479346881e+00, -1.538561538746e+00, -1.700196461237e+00, -1.959316976770e+00, -3.427528679685e+00, -5.272922243974e+00, -6.593995750032e+00, -6.044273233542e+00, -7.629262977162e+00, -1.015069002076e+01, -1.726615367433e+00, -1.135004804706e+00, -2.101067380626e+00, -2.920976388448e+00, -4.379158448292e+00, -2.897738319395e-01, -1.204043168413e+00, -1.567907843610e+00, -1.654987836987e+00, -2.747145752722e+00, -2.744999979469e+00, -8.141849364972e+00, -9.888588149655e+00, -1.017211778784e+01, -1.471319648585e+01, -1.581335685234e+01, -2.655417882343e+01, -2.795297112478e+01, -4.606812517919e+01, -5.251168242716e+01, -5.747281968730e+01, -5.317605197312e+01, -6.434339762912e+01, -7.359481413824e+00, -8.164164536907e+00, -7.926690221317e+00, -9.101347286718e+00, -1.111871433529e+01, -1.671406865188e+01, -1.689231735388e+02, -5.085692898996e-01, -1.870361267917e+00, -3.024155250533e+00, -9.248404616030e+00, -1.353543685892e+01, -1.642530844963e+01, -1.732440024446e+01, -2.309102833890e+01, -2.862502394262e+01, -5.389470429091e-01, -1.401404093684e+00, -5.083936381221e+00, -1.093796113621e+00, -5.018417154523e+00, -1.672134837867e+00, -1.126175538223e+00, -3.461342987826e+00, -7.266035608809e+00, -1.665876292123e+00, -1.287287315941e+00, -2.458367914823e+00, -4.155503207580e+00, -7.262582594804e+00, -1.080175286095e+01, -2.148273722774e+01, -2.118918737437e+01, -2.442766831647e+01, -2.934148853651e-01, -1.198951137691e+00, -5.330540185911e+00, -2.470526290020e+01, -8.331846917000e+01, -8.541310065975e+01, -9.037412006536e+01, -1.504303589770e+02, -2.127773259798e+02, -2.761143201494e+02, -3.396483905922e+02, -8.919174229866e+01, -9.599823708615e+01, -1.038460408955e+02, -1.695860581939e+02, -2.476589981286e+02, -8.818862065268e+01, -9.429363688749e+01, -1.024094645842e+02, -1.190070524906e+02, -3.108575290883e+01, -3.252042819881e+01, -3.899811131342e+01, -6.661083224146e+01, -5.299736762701e+01, -1.861281504789e+00, -1.829869946057e+00, -8.137656316599e+00, -5.223606658528e+01, -3.751643038904e+01, -1.463042322839e+00, -2.031024301520e+00, -2.951222082603e+00, -2.319844437544e+00, -3.891087019464e+00, -6.046893711397e+00, -6.085727442190e+00, -5.950523127845e+00, -1.439379366855e+01, -7.973034132755e+00, -8.011861860246e+00, -8.062543362831e+00, -2.434811423231e+00, -5.506849494178e+00 }};

    Vec<TestDataActivityModelPhreeqc> datalist = { data1, data2, data3 };

    for(auto const& data : datalist)
    {
        // Createt he database object for given database name
        PhreeqcDatabase db(data.dbname);

        // The list of aqueous species in the aqueous solution.
        auto const species = SpeciesList(data.species);

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPhreeqc(db)(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // The temperature used in the PHREEQC script (in K)
        auto const T = data.TdegC + 273.15;

        // The pressure used in the PHREEQC script (in Pa)
        auto const P = data.Patm * 101325.0;

        // Compute the mole fractions of the species
        auto const x = data.n/data.n.sum();

        // Evaluate the activity props function
        fn(props, {T, P, x});

        // Check evaluated activities against the expected values reported by PHREEQC.
        for(auto const& [i, sp] : enumerate(species))
        {
            // Skip species whose original amounts were zero, as their
            // activities reported by PHREEQC are not zero, but most likely
            // resulting from mass-action equation calculations.
            if(data.n[i] == 1e-99)
                continue;

            auto const la_actual = props.ln_a[i].val()/ln10;
            auto const la_expected = data.la[i].val();
            auto const lg_actual = la_actual - log10(data.n[i]);
            auto const lg_expected = la_expected - log10(data.n[i]);

            INFO("Database: " << data.dbname << ", Species: " << sp.name() << ", la_actual: " << la_actual << ", la_expected: " << la_expected << ", lg_actual: " << lg_actual << ", lg_expected: " << lg_expected);
            CHECK( la_actual == Approx(la_expected).epsilon(1e-2) ); // Check Reaktoro's calculated activities (log base 10) are within 1% of PHREEQC's reported values.
        }
    }

    CHECK_THROWS( ActivityModelPhreeqc(PhreeqcDatabase()) ); // Not initialized
    CHECK_THROWS( ActivityModelPhreeqc(PhreeqcDatabase("pitzer.dat")) ); // Pitzer not supported
    CHECK_THROWS( ActivityModelPhreeqc(PhreeqcDatabase("sit.dat")) ); // SIT not supported

    WHEN("The solution is too saline beyond limit point supported by PHREEQC's activity model")
    {
        PhreeqcDatabase db("phreeqc.dat");

        // The list of aqueous species in the aqueous solution.
        auto const species = SpeciesList("H2O H+ OH- Na+ Cl-");

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPhreeqc(db)(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // The temperature used in the PHREEQC script (in K)
        auto const T = 25.0 + 273.15;

        // The pressure used in the PHREEQC script (in Pa)
        auto const P = 1.0 * 101325.0;

        // Compute the mole fractions of the species
        auto const x_with_xw_above_threshold = ArrayXr{{ 0.80, 1e-7, 1e-7, 0.10, 0.10 }}; // xw >= 0.48549234635
        auto const x_with_xw_below_threshold = ArrayXr{{ 0.48, 1e-7, 1e-7, 0.26, 0.26 }}; // xw < 0.48549234635

        // Disable warning about solution being too saline
        Warnings::disable(548);

        // Evaluate the activity props function for both cases of xw
        CHECK_NOTHROW( fn(props, {T, P, x_with_xw_above_threshold}) );
        CHECK_THROWS( fn(props, {T, P, x_with_xw_below_threshold}) );
    }
}
