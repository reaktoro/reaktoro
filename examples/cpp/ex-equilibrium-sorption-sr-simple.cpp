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
//   • Svetlana Kyas (23 November 2021)
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
    AqueousPhase aqueous_phase(speciate("H O Sr"));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOSr+ Hfo_wOSrOH";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g")
               .setMass(4.45, "g");

//    SURFACE 1
//    Hfo_w  1e-3  60   4.45
//    Hfo_s  2.5e-5
//
    // Defined strong site of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w")
               .setAmount(1e-3, "mol");

    // Defined weak site of the complexation surface
    ComplexationSurfaceSite site_Hfo_s;
    site_Hfo_s.setName("Hfo_s")
              .setAmount(0.025e-3, "mol");

    // Add weak site
    surface_Hfo.addSite(site_Hfo_s);

    // Add sorption species
    surface_Hfo.addSurfaceSpecies(list);

    // Add specified surface as parameters for the activity model for the complexation surface
    ActivityModelSurfaceComplexationParams params;
    params.surface = surface_Hfo;

    // Define surface complexation phase and set an activity model
    SurfaceComplexationPhase complexation_phase_Hfo(list_str);
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));

    std::cout << surface_Hfo << std::endl;

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(25.0, "celsius");
    solutionstate.pressure(1.0, "bar");
    solutionstate.set("H2O"    , 1.00, "kg");
    solutionstate.set("Sr+2"  , 1e-6, "mmol");
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

    ComplexationSurfaceProps surface_props(surface_Hfo, solutionstate);
    std::cout << surface_props << std::endl;

    // ------------------------------Surface composition------------------------------
    //
    //Diffuse Double Layer Surface-Complexation Model
    //
    //Hfo
    //	  1.291e-07  Surface charge, eq
    //	  4.664e-05  sigma, C/m�
    //	  5.017e-02  psi, V
    //	 -1.953e+00  -F*psi/RT
    //	  1.419e-01  exp(-F*psi/RT)
    //	  6.000e+01  specific area, m�/g
    //	  2.670e+02  m� for   4.450e+00 g
    //
    //
    //Hfo_s
    //	  2.500e-05  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           1.919e-05       0.768   1.919e-05      -4.717
    //	Hfo_sOH2+         2.906e-06       0.116   2.906e-06      -5.537
    //	Hfo_sO-           2.903e-06       0.116   2.903e-06      -5.537
    //	Hfo_sOHSr+2       3.794e-11       0.000   3.794e-11     -10.421
    //
    //Hfo_w
    //	  1.000e-03  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_wOH           7.676e-04       0.768   7.676e-04      -3.115
    //	Hfo_wOH2+         1.162e-04       0.116   1.163e-04      -3.935
    //	Hfo_wO-           1.161e-04       0.116   1.161e-04      -3.935
    //	Hfo_wOSr+         5.023e-13       0.000   5.023e-13     -12.299
    //	Hfo_wOSrOH        6.176e-16       0.000   6.176e-16     -15.209
    //
    //-----------------------------Solution composition------------------------------
    //
    //	Elements           Molality       Moles
    //
    //	Sr                9.616e-10   9.616e-10

    return 0;
}