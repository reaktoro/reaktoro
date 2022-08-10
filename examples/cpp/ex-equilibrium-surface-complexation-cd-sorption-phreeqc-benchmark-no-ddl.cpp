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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Svetlana Kyas (10 August 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Core/Utils.hpp>

using namespace Reaktoro;

int main()
{
    // Initialize the database
    auto dbphreeqc = PhreeqcDatabase("phreeqc.dat");

    // Define an aqueous phase
    AqueousPhase aqueous_phase(speciate("H O Cd"));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCd+";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

//    String list_str_s = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_sOCd+";
//    String list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCd+";
//    SpeciesList list_s = slist.withNames(StringList(list_str_s));
//    SpeciesList list_w = slist.withNames(StringList(list_str_w));

//    SURFACE 		1
//    Hfo_w 1.0 100 1 # 1.0mol weak site, 100m2/g ssa, 1g Hfo
//    Hfo_s 1.0 # 1.0mol strong site
//    #-no_edl # no ddl

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(100, "m2/g")
               .setMass(1, "g");

    // Defined the sites of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w").setAmount(1.0, "mol");
    surface_Hfo.addSite("Hfo_s", "_s").setAmount(1.0, "mol");

    surface_Hfo.addSurfaceSpecies(list);

    // Add specified surface as parameters for the activity model for the complexation surface
    ActivityModelSurfaceComplexationParams params;
    params.surface = surface_Hfo;

    // Define surface complexation phase and set an activity model
    SurfaceComplexationPhase complexation_phase_Hfo(list_str);
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));

//    SurfaceComplexationPhase complexation_phase_Hfo_w(list_str_w);
//    SurfaceComplexationPhase complexation_phase_Hfo_s(list_str_s);
//    //complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params));
//    complexation_phase_Hfo_w.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));
//    complexation_phase_Hfo_s.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);
    //ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo_w, complexation_phase_Hfo_s);

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(25.0, "celsius");
    solutionstate.pressure(1.0, "atm");
    solutionstate.set("H2O" , 1.0, "kg");
    solutionstate.set("Cd+2", 1.0, "mmol");
    solutionstate.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    solutionstate.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    auto res = solver.solve(solutionstate);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "succeed       = " << res.optima.succeeded << std::endl;
    std::cout << "solutionstate = \n" << solutionstate << std::endl;

    AqueousProps aprops(solutionstate);
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "I  = " << aprops.ionicStrength() << std::endl;

    ChemicalProps props(solutionstate);
    std::cout << "Dissolved(Cd)  = " << props.elementAmountInPhase("Cd", "AqueousPhase") << std::endl;

    ComplexationSurfaceProps surface_props(surface_Hfo, solutionstate);
    std::cout << surface_props << std::endl;

    return 0;

    // -----------------------------Solution composition------------------------------
    //
    //	Elements           Molality       Moles
    //
    //	Cd                1.000e-03   1.000e-03
    //
    //----------------------------Description of solution----------------------------
    //
    //                                       pH  =   7.000
    //                                       pe  =   4.000
    //      Specific Conductance (�S/cm,  25�C)  = 100
    //                          Density (g/cm�)  =   0.99717
    //                               Volume (L)  =   1.00295
    //                        Activity of water  =   1.000
    //                 Ionic strength (mol/kgw)  =   1.999e-03
    //                       Mass of water (kg)  =   1.000e+00
    //                 Total alkalinity (eq/kg)  =   7.230e-07
    //                         Temperature (�C)  =  25.00
    //                  Electrical balance (eq)  =   1.999e-03
    // Percent error, 100*(Cat-|An|)/(Cat+|An|)  =  99.99
    //                               Iterations  =   4
    //                                  Total H  = 1.110124e+02
    //                                  Total O  = 5.550622e+01
    //
    //----------------------------Distribution of species----------------------------
    //
    //                                               Log       Log       Log    mole V
    //   Species          Molality    Activity  Molality  Activity     Gamma    cm�/mol
    //
    //   OH-             1.064e-07   1.012e-07    -6.973    -6.995    -0.022     -4.09
    //   H+              1.047e-07   1.000e-07    -6.980    -7.000    -0.020      0.00
    //   H2O             5.551e+01   1.000e+00     1.744    -0.000     0.000     18.07
    //Cd            1.000e-03
    //   Cd+2            9.993e-04   8.196e-04    -3.000    -3.086    -0.086    -18.71
    //   CdOH+           7.163e-07   6.817e-07    -6.145    -6.166    -0.022     (0)
    //   Cd2OH+3         4.274e-09   2.736e-09    -8.369    -8.563    -0.194     (0)
    //   Cd(OH)2         3.659e-10   3.661e-10    -9.437    -9.436     0.000     (0)
    //   Cd(OH)3-        4.316e-16   4.107e-16   -15.365   -15.386    -0.022     (0)
    //   Cd(OH)4-2       4.463e-23   3.661e-23   -22.350   -22.436    -0.086     (0)
    //H(0)          1.415e-25
    //   H2              7.076e-26   7.079e-26   -25.150   -25.150     0.000     28.61
    //O(0)          0.000e+00
    //   O2              0.000e+00   0.000e+00   -42.080   -42.080     0.000     30.40
    //
    //------------------------------Saturation indices-------------------------------
    //
    //  Phase               SI** log IAP   log K(298 K,   1 atm)
    //
    //  Cd(OH)2          -2.74     10.91   13.65  Cd(OH)2
    //  H2(g)           -22.05    -25.15   -3.10  H2
    //  H2O(g)           -1.50     -0.00    1.50  H2O
    //  O2(g)           -39.19    -42.08   -2.89  O2
    //
    //**For a gas, SI = log10(fugacity). Fugacity = pressure * phi / 1 atm.
    //  For ideal gases, phi = 1.
    //
    //-----------------------------------------
    //Beginning of batch-reaction calculations.
    //-----------------------------------------
    //
    //Reaction step 1.
    //
    //Using solution 1.
    //Using surface 1.
    //
    //----------------------------------User print-----------------------------------
    //
    //Sorbed(Cd) =    1.0374e-03
    //Total(Cd)  =    3.4852e-12
    //
    //------------------------------Surface composition------------------------------
    //
    //Hfo
    //	  2.001e-03  Surface charge, eq
    //Hfo_s
    //	  1.000e+00  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           7.669e-01       0.767   7.955e-01      -0.099
    //	Hfo_sOH2+         1.163e-01       0.116   1.207e-01      -0.918
    //	Hfo_sO-           1.158e-01       0.116   1.201e-01      -0.920
    //	Hfo_sOCd+         9.996e-04       0.001   1.037e-03      -2.984
    //
    //Hfo
    //	  2.001e-03  Surface charge, eq
    //Hfo_w
    //	  1.000e+00  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_wOH           7.676e-01       0.768   7.963e-01      -0.099
    //	Hfo_wOH2+         1.164e-01       0.116   1.208e-01      -0.918
    //	Hfo_wO-           1.159e-01       0.116   1.203e-01      -0.920
    //	Hfo_wOCd+         4.171e-07       0.000   4.327e-07      -6.364
    //
    //-----------------------------Solution composition------------------------------
    //
    //	Elements           Molality       Moles
    //
    //	Cd                3.485e-12   3.360e-12
    //
    //----------------------------Description of solution----------------------------
    //
    //                                       pH  =   8.109      Charge balance
    //                                       pe  =  -9.692      Adjusted to redox equilibrium
    //      Specific Conductance (�S/cm,  25�C)  = 0
    //                          Density (g/cm�)  =   0.97042
    //                               Volume (L)  =   0.99543
    //                        Activity of water  =   0.982
    //                 Ionic strength (mol/kgw)  =   6.436e-07
    //                       Mass of water (kg)  =   9.640e-01
    //                 Total alkalinity (eq/kg)  =   1.272e-06
    //                         Temperature (�C)  =  25.00
    //                  Electrical balance (eq)  =  -1.226e-06
    // Percent error, 100*(Cat-|An|)/(Cat+|An|)  = -98.79
    //                               Iterations  =  20
    //                                  Total H  = 1.090124e+02
    //                                  Total O  = 5.350622e+01
    //
    //----------------------------Distribution of species----------------------------
    //
    //                                               Log       Log       Log    mole V
    //   Species          Molality    Activity  Molality  Activity     Gamma    cm�/mol
    //
    //   OH-             1.279e-06   1.278e-06    -5.893    -5.893    -0.000     -4.14
    //   H+              7.787e-09   7.779e-09    -8.109    -8.109    -0.000      0.00
    //   H2O             5.551e+01   9.824e-01     1.744    -0.008     0.000     18.07
    //Cd            3.485e-12
    //   Cd+2            3.449e-12   3.436e-12   -11.462   -11.464    -0.002    -18.86
    //   CdOH+           3.612e-14   3.609e-14   -13.442   -13.443    -0.000     (0)
    //   Cd(OH)2         2.447e-16   2.447e-16   -15.611   -15.611     0.000     (0)
    //   Cd(OH)3-        3.471e-21   3.468e-21   -20.460   -20.460    -0.000     (0)
    //   Cd2OH+3         6.125e-25   6.073e-25   -24.213   -24.217    -0.004     (0)
    //   Cd(OH)4-2       3.918e-27   3.903e-27   -26.407   -26.409    -0.002     (0)
    //H(0)          2.075e+00
    //   H2              1.037e+00   1.037e+00     0.016     0.016     0.000     28.61
    //O(0)          0.000e+00
    //   O2              0.000e+00   0.000e+00   -92.427   -92.427     0.000     30.40
    //
    //------------------------------Saturation indices-------------------------------
    //
    //  Phase               SI** log IAP   log K(298 K,   1 atm)
    //
    //  Cd(OH)2          -8.91      4.74   13.65  Cd(OH)2
    //  H2(g)             3.12      0.02   -3.10  H2
    //  H2O(g)           -1.51     -0.01    1.50  H2O
    //  O2(g)           -89.53    -92.43   -2.89  O2
}