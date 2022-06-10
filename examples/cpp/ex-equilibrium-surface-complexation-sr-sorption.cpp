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
    AqueousPhase aqueous_phase(speciate("H O Cl Ca Sr"));
    aqueous_phase.setActivityModel(ActivityModelHKF());


//    Hfo_sOH           1.126e-05       0.450   1.126e-05      -4.949
//    Hfo_sOHCa+2       1.034e-05       0.413   1.034e-05      -4.986
//    Hfo_sOH2+         1.714e-06       0.069   1.714e-06      -5.766
//    Hfo_sO-           1.693e-06       0.068   1.693e-06      -5.771
//    Hfo_sOHSr+2       1.134e-11       0.000   1.134e-11     -10.945
//
//    Hfo_wOH           7.666e-04       0.767   7.666e-04      -3.115
//    Hfo_wOH2+         1.167e-04       0.117   1.167e-04      -3.933
//    Hfo_wO-           1.153e-04       0.115   1.153e-04      -3.938
//    Hfo_wOCa+         1.364e-06       0.001   1.364e-06      -5.865
//    Hfo_wOSr+         2.542e-13       0.000   2.542e-13     -12.595
//    Hfo_wOSrOH        3.108e-16       0.000   3.108e-16     -15.508


    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";
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
    surface_Hfo.addSite(site_Hfo_s);

    surface_Hfo.addSurfaceSpecies(list);

    std::cout << surface_Hfo << std::endl;

    // Add specified surface as parameters for the activity model for the complexation surface
    ActivityModelSurfaceComplexationParams params;
    params.surface = surface_Hfo;

    // Define surface complexation phase and set an activity model
    SurfaceComplexationPhase complexation_phase_Hfo(list_str);
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in bar

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(T, "celsius");
    solutionstate.pressure(P, "bar");
    solutionstate.set("H2O"    , 1.00, "kg");
//    Ca     1
//    Cl     2
//    Sr     1e-006
    solutionstate.set("Cl-"  , 2e+0, "mmol");
    solutionstate.set("Ca+2"  , 1e+0, "mmol");
    solutionstate.set("Sr+2"  , 1e-6, "mmol");
//    Hfo_w 1e-3 60 4.45 # 1e-3mol weak site, 60m2/g s.spec, 4.45g ferrihyd
//    Hfo_s 0.025e-3 # 0.025e-3mol strong site
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

    return 0;
}