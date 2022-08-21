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
    auto elements = "H O Na Cl";
    AqueousPhase aqueous_phase(speciate(elements));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define ion exchange species list
    String list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- "
                      "Hfo_wOH Hfo_wOH2+ Hfo_wO-";
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    SpeciesList list = slist.withNames(StringList(list_str));

    String list_str_s = "Hfo_sOH Hfo_sOH2+ Hfo_sO-";
    String list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO-";
//    SURFACE 		1
//    Hfo_w 1.0 100 1 # 1.0mol weak site, 100m2/g s.spec, 10g ferrihyd
//    Hfo_s 1.0 # 1.0mol strong site
//    #-no_edl # no ddl

    // Create complexation surface
//    ComplexationSurface surface_Hfo("Hfo");
//    surface_Hfo.setSpecificSurfaceArea(100, "m2/g")
//               .setMass(1, "g");

    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g");

    // Defined the sites of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol");
    surface_Hfo.addSite("Hfo_s", "_s").setAmount(0.025e-3, "mol");

    // Add species to the surface and corresponding sites
    surface_Hfo.addSurfaceSpecies(list);

    // Add specified surface as parameters for the activity model for the complexation surface
    ActivityModelSurfaceComplexationSiteParams params_site;
    params_site.surface = surface_Hfo;

    // Define surface complexation phase and set an activity model
    params_site.site_tag = "_w";
    SurfaceComplexationPhase complexation_phase_Hfo_w(list_str_w);
    complexation_phase_Hfo_w.setName("Hfo_w");
    complexation_phase_Hfo_w.setActivityModel(ActivityModelSurfaceComplexationSiteNoDDL(params_site));

    params_site.site_tag = "_s";
    SurfaceComplexationPhase complexation_phase_Hfo_s(list_str_s);
    complexation_phase_Hfo_s.setName("Hfo_s");
    complexation_phase_Hfo_s.setActivityModel(ActivityModelSurfaceComplexationSiteNoDDL(params_site));

    ActivityModelDDLParams params_dll;
    params_dll.output = true;

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, complexation_phase_Hfo_w, complexation_phase_Hfo_s);

    auto pHs = {4, 5, 6, 7, 8, 9, 10};
    //auto pHs = {4};

    ChemicalProps props(system);
    ComplexationSurfaceSiteProps site_w_props(surface_Hfo.sites()["_w"], system);
    ComplexationSurfaceSiteProps site_s_props(surface_Hfo.sites()["_s"], system);

    EquilibriumOptions opts;
    //opts.optima.output.active = true;

    for(auto pH : pHs) {
        // Define initial equilibrium state
        ChemicalState solutionstate(system);
        solutionstate.set("H2O", 1.0, "kg");
        solutionstate.set("Na+", 1e-1, "mol");
        solutionstate.set("Cl-", 1e-1, "mol");
        solutionstate.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol");
        solutionstate.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol");

        EquilibriumSpecs specs(system);
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumConditions conditions(specs);
        conditions.temperature(25.0, "celsius");
        conditions.pressure(1.0, "atm");
        conditions.pH(pH);

        // Define equilibrium solver and corresponding equilibrium options
        EquilibriumSolver solver(specs);
        solver.setOptions(opts);

        // Equilibrate given initial state
        auto res = solver.solve(solutionstate, conditions);
//        std::cout << "*******************************************" << std::endl;
//        std::cout << "After equilibration: " << std::endl;
//        std::cout << "*******************************************" << std::endl;
//        std::cout << "succeed       = " << res.optima.succeeded << std::endl;
//        std::cout << "solutionstate = \n" << solutionstate << std::endl;

//        AqueousProps aprops(solutionstate);
//        std::cout << "pH = " << aprops.pH() << std::endl;
//        std::cout << "I  = " << aprops.ionicStrength() << std::endl;
//        std::cout << "aprops = \n" << aprops << std::endl;

        props.update(solutionstate);

        //ComplexationSurfaceProps surface_props(surface_Hfo, solutionstate);
        site_w_props.update(solutionstate);
        site_s_props.update(solutionstate);

        //std::cout << "surface_props = \n" << surface_props << std::endl;

        std::cout << "pH = " << pH
                  << ", x(Hfo_wOH) = " << site_w_props.speciesFraction("Hfo_wOH")
                  << ", x(Hfo_wOH2+) = " << site_w_props.speciesFraction("Hfo_wOH2+")
                  << ", x(Hfo_wO-) = " << site_w_props.speciesFraction("Hfo_wO-")
                  << ", x(Hfo_sOH) = " << site_s_props.speciesFraction("Hfo_sOH")
                  << ", x(Hfo_sOH2+) = " << site_s_props.speciesFraction("Hfo_sOH2+")
                  << ", x(Hfo_sO-) = " << site_s_props.speciesFraction("Hfo_sO-") << std::endl;
    }
    return 0;

}