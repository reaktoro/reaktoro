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
    AqueousPhase aqueous_phase(speciate("H O Cd Cl"));
    aqueous_phase.setActivityModel(chain(ActivityModelHKF(), ActivityModelElectrostatics(ActivityModelDDLParams())));

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCd+";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

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
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params));

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(25.0, "celsius");
    solutionstate.pressure(1.0, "atm");
    solutionstate.set("H2O" , 1.0, "kg");
    solutionstate.set("Cd+2", 1.0, "mmol");
    solutionstate.set("Cl-" , 2.0, "mmol"); // + 1e-6 is to balance the charge of Sr+2
    solutionstate.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    solutionstate.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    EquilibriumOptions opts;
//    opts.optima.output.active = true;
//    opts.optima.maxiters = 300;

    // Define equilibrium solver
    EquilibriumSolver solver(system);
    solver.setOptions(opts);

    // Equilibrate given initial state
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
    std::cout << "dissolved(Cd)   = " << props.elementAmountInPhase("Cd", "AqueousPhase") << std::endl;
    std::cout << "charge(aqueous) = " <<  props.chargeInPhase("AqueousPhase") << std::endl;

    ComplexationSurfaceProps surface_props(surface_Hfo, solutionstate);
    std::cout << "surface_props = \n" << surface_props << std::endl;

    return 0;

    // ----------------------------------User print-----------------------------------
    //
    //mol(Cd+)   =    7.3788e-04
    //Sorbed(Cd) =    2.9949e-04
    //Total(Cd)  =    7.3788e-04
    //
    //------------------------------Surface composition------------------------------
    //
    //Diffuse Double Layer Surface-Complexation Model
    //
    //Hfo
    //	  3.816e-04  Surface charge, eq
    //	  3.683e-01  sigma, C/m�
    //	  2.602e-01  psi, V
    //	 -1.013e+01  -F*psi/RT
    //	  4.001e-05  exp(-F*psi/RT)
    //	  1.000e+02  specific area, m�/g
    //	  1.000e+02  m� for   1.000e+00 g
    //
    //
    //Hfo_s
    //	  1.000e+00  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           7.674e-01       0.767   7.961e-01      -0.099
    //	Hfo_sOH2+         1.162e-01       0.116   1.205e-01      -0.919
    //	Hfo_sO-           1.161e-01       0.116   1.205e-01      -0.919
    //	Hfo_sOCd+         2.886e-04       0.000   2.994e-04      -3.524
    //
    //Hfo_w
    //	  1.000e+00  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_wOH           7.676e-01       0.768   7.963e-01      -0.099
    //	Hfo_wOH2+         1.162e-01       0.116   1.206e-01      -0.919
    //	Hfo_wO-           1.162e-01       0.116   1.205e-01      -0.919
    //	Hfo_wOCd+         1.203e-07       0.000   1.248e-07      -6.904
    //
    //-----------------------------Solution composition------------------------------
    //
    //	Elements           Molality       Moles
    //
    //	Cd                7.379e-04   7.113e-04
    //
    //----------------------------Description of solution----------------------------
    //
    //                                       pH  =   3.712      Charge balance
    //                                       pe  =  -5.295      Adjusted to redox equilibrium
    //      Specific Conductance (�S/cm,  25�C)  = 144
    //                          Density (g/cm�)  =   0.97051
    //                               Volume (L)  =   0.99542
    //                        Activity of water  =   0.982
    //                 Ionic strength (mol/kgw)  =   1.577e-03
    //                       Mass of water (kg)  =   9.640e-01
    //                 Total alkalinity (eq/kg)  =  -2.023e-04
    //                         Temperature (�C)  =  25.00
    //                  Electrical balance (eq)  =   1.618e-03
    // Percent error, 100*(Cat-|An|)/(Cat+|An|)  = 100.00
    //                               Iterations  =  24
    //                                  Total H  = 1.090126e+02
    //                                  Total O  = 5.350622e+01
    //
    //----------------------------Distribution of species----------------------------
    //
    //                                               Log       Log       Log    mole V
    //   Species          Molality    Activity  Molality  Activity     Gamma    cm�/mol
    //
    //   H+              2.023e-04   1.941e-04    -3.694    -3.712    -0.018      0.00
    //   OH-             5.357e-11   5.124e-11   -10.271   -10.290    -0.019     -4.10
    //   H2O             5.551e+01   9.823e-01     1.744    -0.008     0.000     18.07
    //Cd            7.379e-04
    //   Cd+2            7.379e-04   6.181e-04    -3.132    -3.209    -0.077    -18.72
    //   CdOH+           2.720e-10   2.602e-10    -9.565    -9.585    -0.019     (0)
    //   Cd2OH+3         1.174e-12   7.878e-13   -11.930   -12.104    -0.173     (0)
    //   Cd(OH)2         7.072e-17   7.074e-17   -16.150   -16.150     0.000     (0)
    //   Cd(OH)3-        4.200e-26   4.018e-26   -25.377   -25.396    -0.019     (0)
    //   Cd(OH)4-2       2.164e-36   1.813e-36   -35.665   -35.742    -0.077     (0)
    //H(0)          2.075e+00
    //   H2              1.037e+00   1.038e+00     0.016     0.016     0.000     28.61
    //O(0)          0.000e+00
    //   O2              0.000e+00   0.000e+00   -92.428   -92.428     0.000     30.40
    //
    //------------------------------Saturation indices-------------------------------
    //
    //  Phase               SI** log IAP   log K(298 K,   1 atm)
    //
    //  Cd(OH)2          -9.45      4.20   13.65  Cd(OH)2
    //  H2(g)             3.12      0.02   -3.10  H2
    //  H2O(g)           -1.51     -0.01    1.50  H2O
    //  O2(g)           -89.54    -92.43   -2.89  O2
}