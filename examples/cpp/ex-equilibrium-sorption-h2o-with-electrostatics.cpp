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
    auto elements = "H O";
    AqueousPhase aqueous_phase(speciate(elements));
    aqueous_phase.setActivityModel(chain(ActivityModelHKF(), ActivityModelElectrostatics(ActivityModelDDLParams())));

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO-";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

//    SURFACE 		1
//    Hfo_w 1.0 100 1 # 1.0mol weak site, 100m2/g s.spec, 10g ferrihyd
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
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params));

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo);

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.temperature(25.0, "celsius");
    solutionstate.pressure(1.0, "atm");
    solutionstate.set("H2O" , 1.0, "kg");
    solutionstate.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    solutionstate.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");
//    auto scale = 1e-6;
//    solutionstate.set("H2O!" , scale*1.0, "kg");
//    solutionstate.set("Cl-!" , scale*2e+0, "mmol"); // + 1e-6 is to balance the charge of Sr+2
//    solutionstate.set("Cd+2!", scale*1.0, "mmol");

//    std::cout << "*******************************************" << std::endl;
//    std::cout << "Before equilibration: " << std::endl;
//    std::cout << "*******************************************" << std::endl;
    ChemicalProps props(solutionstate);
    //std::cout << "aq.phase charge = " << props.chargeInPhase("AqueousPhase") << std::endl;
    //std::cout << "props = \n" << props << std::endl;

    EquilibriumOptions opts;
    opts.optima.output.active = true;
//    opts.optima.maxiters = 250;

    // Define equilibrium solver and corresponding equilibrium options
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
    props.update(solutionstate);
    std::cout << "aq.phase charge = " <<  props.chargeInPhase("AqueousPhase") << std::endl;
    std::cout << "pH = " << aprops.pH() << std::endl;
    std::cout << "I  = " << aprops.ionicStrength() << std::endl;

    std::cout << "a(H+)  = " << props.speciesActivity("H+") << std::endl;
    std::cout << "a(OH-)  = " << props.speciesActivity("OH-") << std::endl;
    std::cout << "gamma(H+)  = " << props.speciesActivityCoefficient("H+") << std::endl;
    std::cout << "gamma(OH-)  = " << props.speciesActivityCoefficient("OH-") << std::endl;

    std::cout << "aprops = \n" << aprops << std::endl;

    ComplexationSurfaceProps surface_props(surface_Hfo, solutionstate);
    std::cout << "surface_props = \n" << surface_props << std::endl;
    std::cout << "succeed       = " << res.optima.succeeded << std::endl;

    return 0;
}