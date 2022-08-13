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
//   • Svetlana Kyas (15 June 2022)
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
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+ ";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g")
               .setMass(4.45, "g");

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

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(25.0, "celsius");
    conditions.pressure(1.0, "bar");
    conditions.pH(6.0);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.set("H2O"   , 1.00, "kg");
    state.set("Cl-"   , 2e+0, "mmol");
    state.set("Ca+2"  , 1e+0, "mmol");
    state.set("Cd+2"  , 1e-6, "mmol");
    state.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    state.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    // Define equilibrium solver and equilibrate given initial state
    auto res = solver.solve(state, conditions);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "succeed       = " << res.optima.succeeded << std::endl;
    std::cout << "solutionstate = \n" << state << std::endl;

    AqueousProps aprops(state);
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "aprops = \n" << aprops << std::endl;

    ComplexationSurfaceProps surface_props(surface_Hfo, state);
    std::cout << surface_props << std::endl;

    ChemicalProps props(state);
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "Kd = " << (state.speciesAmount("Hfo_sOCd+") + state.speciesAmount("Hfo_wOCd+")) / props.elementAmount("Cd") << std::endl;

    // Diffuse Double Layer Surface-Complexation Model
    //
    //Hfo
    //	  2.348e-05  Surface charge, eq
    //	  8.485e-03  sigma, C/m�
    //	  5.621e-02  psi, V
    //	 -2.188e+00  -F*psi/RT
    //	  1.122e-01  exp(-F*psi/RT)
    //	  6.000e+01  specific area, m�/g
    //	  2.670e+02  m� for   4.450e+00 g
    //
    //
    //Hfo_s
    //	  2.500e-05  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           1.126e-05       0.450   1.126e-05      -4.949
    //	Hfo_sOHCa+2       1.033e-05       0.413   1.033e-05      -4.986
    //	Hfo_sOH2+         1.714e-06       0.069   1.714e-06      -5.766
    //	Hfo_sO-           1.693e-06       0.068   1.693e-06      -5.771
    //	Hfo_sOCd+         9.472e-10       0.000   9.472e-10      -9.024
    //
    //Hfo_w
    //	  1.000e-03  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_wOH           7.666e-04       0.767   7.666e-04      -3.115
    //	Hfo_wOH2+         1.167e-04       0.117   1.167e-04      -3.933
    //	Hfo_wO-           1.153e-04       0.115   1.153e-04      -3.938
    //	Hfo_wOCa+         1.364e-06       0.001   1.364e-06      -5.865
    //	Hfo_wOCd+         2.689e-11       0.000   2.689e-11     -10.570
    //
    //-----------------------------Solution composition------------------------------
    //
    //	Elements           Molality       Moles
    //
    //	Ca                9.883e-04   9.883e-04
    //	Cd                2.593e-11   2.593e-11
    //	Cl                2.000e-03   2.000e-03
    // -----------------------------------------------------------------------------
    //
    // Hfo
    //	  8.384e-05  Surface charge, eq
    //	  3.030e-02  sigma, C/m�
    //	  1.159e-01  psi, V
    //	 -4.510e+00  -F*psi/RT
    //	  1.100e-02  exp(-F*psi/RT)
    //	  6.000e+01  specific area, m�/g
    //	  2.670e+02  m� for   4.450e+00 g
    //
    //
    //Hfo_s
    //	  2.500e-05  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           1.879e-05       0.752   1.879e-05      -4.726
    //	Hfo_sOH2+         4.031e-06       0.161   4.031e-06      -5.395
    //	Hfo_sO-           2.007e-06       0.080   2.007e-06      -5.697
    //	Hfo_sOHCa+2       1.678e-07       0.007   1.678e-07      -6.775
    //	Hfo_sOCd+         4.169e-10       0.000   4.169e-10      -9.380
    //
    //Hfo_w
    //	  1.000e-03  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_wOH           7.568e-04       0.757   7.568e-04      -3.121
    //	Hfo_wOH2+         1.623e-04       0.162   1.623e-04      -3.790
    //	Hfo_wO-           8.084e-05       0.081   8.084e-05      -4.092
    //	Hfo_wOCa+         9.297e-09       0.000   9.297e-09      -8.032
    //	Hfo_wOCd+         6.998e-12       0.000   6.998e-12     -11.155
    return 0;
}