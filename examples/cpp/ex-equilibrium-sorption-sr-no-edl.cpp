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
//   • Svetlana Kyas (23 August 2021)
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
    auto elements = "H O Cl Ca Sr";
    AqueousPhase aqueous_phase(speciate(elements));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Fetch all the species with Adsorbed aggregate state
    SpeciesList adsorbed_species = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);

    // Defined and add surface species
    String selected_species_s = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2";
    String selected_species_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";

    // Select the names of considered absorbed species
    SpeciesList species_list_s = adsorbed_species.withNames(selected_species_s);
    SpeciesList species_list_w = adsorbed_species.withNames(selected_species_w);
    SpeciesList species_list = species_list_s + species_list_w;

    // Create the surface
    Surface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g");

    // Add species to the surface
    surface_Hfo.addSurfaceSpecies(species_list);

    // Defined strong site of the surface
    surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol");

    // Defined weak site of the surface
    SurfaceSite site_Hfo_s;
    site_Hfo_s.setName("Hfo_s").setAmount(0.025e-3, "mol");
    surface_Hfo.addSite(site_Hfo_s);

    // Print out the surface structure
    std::cout << surface_Hfo << std::endl;

    // Define surface sites as sites phases
    SurfacePhase hfo_w_phase(selected_species_w);
    hfo_w_phase.setName("Hfo_w");
    SurfacePhase hfo_s_phase(selected_species_s);
    hfo_s_phase.setName("Hfo_s");

    // Define parameters for the activity model of the surface and set corresponding activity model
    ActivityModelSorptionParams params_site;
    params_site.surface = surface_Hfo;
    params_site.site_tag = "_w";
    hfo_w_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site));
    params_site.site_tag = "_s";
    hfo_s_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site));

    // Construct the chemical system
    ChemicalSystem system(dbphreeqc, aqueous_phase, hfo_w_phase, hfo_s_phase);

    // Define temperature and pressure
    const auto T = 25.0;
    const auto P = 1.0;

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumOptions opts;
    //opts.optima.output.active = true;
    EquilibriumSolver solver(specs);
    solver.setOptions(opts);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.set("H2O" , 1.00, "kg");
    state.set("Cl-" , 2e+0, "mmol"); // + 1e-6 is to balance the charge of Sr+2
    state.set("Ca+2", 1e+0, "mmol");
    state.set("Sr+2", 1e-3, "mmol");
    // Set amounts of the surface species in the chemical state
    state.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    state.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    std::cout << "*******************************************" << std::endl;
    std::cout << "Before equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    ChemicalProps props(state);
    std::cout << "Aq.phase charge = " << props.chargeInPhase("AqueousPhase") << std::endl;

    // Define equilibrium solver and equilibrate given initial state
    auto res = solver.solve(state, conditions);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "Convergence succeed: " << res.optima.succeeded << std::endl;
    std::cout << "State \n" << state << std::endl;

    AqueousProps aprops(state);
    props.update(state);
    std::cout << "Aqueous properties:" << std::endl;
    std::cout << aprops << std::endl;
    std::cout << "Aq. phase charge = " << props.chargeInPhase("AqueousPhase") << std::endl;
    std::cout << "Sr sorbed        = " << props.elementAmountInPhase("Sr", "Hfo_s")
                                          + props.elementAmountInPhase("Sr", "Hfo_w") << std::endl;
    std::cout << "Sr dissolved     = " << props.elementAmountInPhase("Sr", "AqueousPhase") << std::endl;

    // Evaluate properties of the surface and its sites
    SurfaceProps sprops(surface_Hfo, state);
    std::cout << sprops << std::endl;

    return 0;
}