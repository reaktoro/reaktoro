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
    AqueousPhase aqueous_phase(speciate("H O Cl Ca Cd"));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g")
               .setMass(4.45, "g");

//    SURFACE 		1
//    #2790 ppm Fe /55.85 = 50 mmol/kg*89 =4.45g "Ferrihyd"/kg
//    Hfo_w 1e-3 60 4.45 # 1e-3mol weak site, 60m2/g s.spec, 4.45g ferrihyd
//    Hfo_s 0.025e-3 # 0.025e-3mol strong site
//    -equil  1

    // Defined weak site of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w")
               .setAmount(1e-3, "mol");

    // Defined weak site of the complexation surface
    ComplexationSurfaceSite site_Hfo_s;
    site_Hfo_s.setName("Hfo_s")
              .setAmount(0.025e-3, "mol");
    surface_Hfo.addSite(site_Hfo_s);

    surface_Hfo.addSurfaceSpecies(list);

    // Add specified surface as parameters for the activity model for the complexation surface
    ActivityModelSurfaceComplexationParams params;
    params.surface = surface_Hfo;

    // Define surface complexation phase and set an activity model
    SurfaceComplexationPhase complexation_phase_Hfo(list_str);
    //complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params));
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationGainesThomas(params));

    std::cout << surface_Hfo << std::endl;

    //    Hfo
//    2.844e-05  Surface charge, eq
//    1.028e-02  sigma, C/m�
//    6.416e-02  psi, V
//   -2.497e+00  -F*psi/RT
//    8.231e-02  exp(-F*psi/RT)
//    6.000e+01  specific area, m�/g
//    2.670e+02  m� for   4.450e+00 g
//
//
//        Hfo_s
//    2.500e-05  moles
//    Mole                     Log
//    Species               Moles    Fraction    Molality    Molality
//
//    Hfo_sOH           1.385e-05       0.554   1.385e-05      -4.858
//    Hfo_sOHCa+2       6.925e-06       0.277   6.925e-06      -5.160
//    Hfo_sOH2+         2.223e-06       0.089   2.223e-06      -5.653
//    Hfo_sO-           1.977e-06       0.079   1.977e-06      -5.704
//    Hfo_sOCd+         2.297e-08       0.001   2.297e-08      -7.639
//
//    Hfo_w
//    1.000e-03  moles
//    Mole                     Log
//    Species               Moles    Fraction    Molality    Molality
//
//    Hfo_wOH           7.668e-04       0.767   7.668e-04      -3.115
//    Hfo_wOH2+         1.231e-04       0.123   1.231e-04      -3.910
//    Hfo_wO-           1.094e-04       0.109   1.094e-04      -3.961
//    Hfo_wOCa+         7.049e-07       0.001   7.049e-07      -6.152
//    Hfo_wOCd+         5.301e-10       0.000   5.301e-10      -9.276

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);

    const auto T = 25.0; // temperature in celsius
    const auto P = 1.0;  // pressure in bar

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.setTemperature(T, "celsius");
    solutionstate.setPressure(P, "bar");
    solutionstate.setSpeciesMass("H2O"    , 1.00, "kg");
//    Ca     1
//    Cd   1e-006
//    Cl      2
    solutionstate.setSpeciesAmount("Cl-"  , 2e+0, "mmol");
    solutionstate.setSpeciesAmount("Ca+2"  , 1e+0, "mmol");
    solutionstate.setSpeciesAmount("Cd+2"  , 1e-6, "mmol");
//    Hfo_w 1e-3 60 4.45 # 1e-3mol weak site, 60m2/g s.spec, 4.45g ferrihyd
//    Hfo_s 0.025e-3 # 0.025e-3mol strong site
    solutionstate.setSpeciesAmount("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    solutionstate.setSpeciesAmount("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

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