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
    AqueousPhase aqueous_phase(speciate("H O Cd"));
    aqueous_phase.setActivityModel(ActivityModelHKF());

    // Define surface species
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);
    String list_str_s = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_sOCd+";
    String list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCd+";
    SpeciesList list_s = slist.withNames(StringList(list_str_s));
    SpeciesList list_w = slist.withNames(StringList(list_str_w));
    SpeciesList list = list_s + list_w;

    // Create surface
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
    state.set("Cd+2", 1.0, "mmol");
    state.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol");
    state.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    auto res = solver.solve(state);
    std::cout << "Convergence succeed: " << res.optima.succeeded << std::endl;
    std::cout << "State: \n" << state << std::endl;

    // Evaluate some aqueous and chemical properties
    AqueousProps aprops(state);
    ChemicalProps props(state);
    std::cout << "Aqueous properties:" << std::endl;
    std::cout << "pH           = " << aprops.pH() << std::endl;
    std::cout << "I            = " << aprops.ionicStrength() << std::endl;
    std::cout << "Cd sorbed    = " << props.elementAmountInPhase("Cd", "Hfo_s") + props.elementAmountInPhase("Cd", "Hfo_w") << std::endl;
    std::cout << "Cd dissolved = " << props.elementAmountInPhase("Cd", "AqueousPhase") << std::endl;

    // Evaluate properties of the surface and its sites
    SurfaceProps sprops(surface_Hfo, state);
    std::cout << sprops << std::endl;

    return 0;
}