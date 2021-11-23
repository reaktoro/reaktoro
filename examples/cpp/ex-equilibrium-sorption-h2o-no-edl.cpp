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
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define ion exchange species list
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    String list_str_s = "Hfo_sOH Hfo_sOH2+ Hfo_sO-";
    String list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO-";
    SpeciesList list_s = slist.withNames(StringList(list_str_s));
    SpeciesList list_w = slist.withNames(StringList(list_str_w));
    SpeciesList list = list_s + list_w;

    // Create the surface
    Surface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(100, "m2/g").setMass(1, "g");

    // Defined the sites of the surface
    surface_Hfo.addSite("Hfo_w", "_w").setAmount(1.0, "mol");
    surface_Hfo.addSite("Hfo_s", "_s").setAmount(1.0, "mol");

    // Add species to the surface
    surface_Hfo.addSurfaceSpecies(list);

    std::cout << surface_Hfo << std::endl;

    // Define surface sites as sites phases
    SurfacePhase hfo_w_phase(list_str_w);
    hfo_w_phase.setName("Hfo_w");
    SurfacePhase hfo_s_phase(list_str_s);
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

    // Define initial equilibrium state
    ChemicalState state(system);
    state.temperature(25.0, "celsius");
    state.pressure(1.0, "atm");
    state.set("H2O" , 1.0, "kg");
    state.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    state.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    // Define chemical properties
    ChemicalProps props(state);

    std::cout << "*******************************************" << std::endl;
    std::cout << "Before equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "Aq.phase charge = " <<  props.chargeInPhase("AqueousPhase") << std::endl;

    EquilibriumOptions opts;
    //  opts.optima.output.active = true;

    // Define equilibrium solver and corresponding equilibrium options
    EquilibriumSolver solver(system);
    solver.setOptions(opts);

    // Equilibrate given initial state
    auto res = solver.solve(state);
    std::cout << "*******************************************" << std::endl;
    std::cout << "After equilibration: " << std::endl;
    std::cout << "*******************************************" << std::endl;
    std::cout << "Convergence succeed: " << res.optima.succeeded << std::endl;
    std::cout << "State: \n" << state << std::endl;

    // Evaluate some aqueous and chemical properties
    AqueousProps aprops(state);
    props.update(state);
    std::cout << "Aqueous properties:" << std::endl;
    std::cout << "pH              = " << aprops.pH() << std::endl;
    std::cout << "I               = " << aprops.ionicStrength() << std::endl;
    std::cout << "Aq.phase charge = " <<  props.chargeInPhase("AqueousPhase") << std::endl;
    std::cout << "Aq.phase mass   = " <<  props.phaseProps("AqueousPhase").mass() << std::endl;

    // Evaluate properties of the surface and its sites
    SurfaceProps sprops(surface_Hfo, state);
    std::cout << sprops << std::endl;

    return 0;
}