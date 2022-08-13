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
//   • Svetlana Kyas (14 June 2022)
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

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";
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
    ChemicalState solutionstate(system);
    solutionstate.set("H2O"    , 1.00, "kg");
    solutionstate.set("Cl-"  , 2e+0, "mmol");
    solutionstate.set("Ca+2"  , 1e+0, "mmol");
    solutionstate.set("Sr+2"  , 1e-6, "mmol");
    solutionstate.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    solutionstate.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    // Define equilibrium solver and equilibrate given initial state
    auto res = solver.solve(solutionstate, conditions);
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