// Reaktoro is a unified framework for modeling chemically reactive systems.
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

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

using namespace Reaktoro;

// Reactive transport test includes
#include <demos/cpp/DemosUtils.h>

auto runKinetics(KineticPathParams &) -> void;

int main()
{
    int minute(60);

    KineticPathParams params = {};

    params.t0 = 0;
    params.tfinal = 30 * minute;

    runKinetics(params);
}

auto runKinetics(KineticPathParams & params) -> void{

    // Create results file
    auto filename = params.makeResultsFile("calcite-dissolution-different-ssarea");

    std::vector<double> surface_areas = {1, 10, 100, 1000};
    std::vector<ChemicalState> chemical_states;
    const unsigned int states_num = surface_areas.size();

    ChemicalEditor editor;
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
    editor.addMineralPhase("Calcite");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");

    for(unsigned i = 0; i < states_num; ++i)
    {
        reaction.setSpecificSurfaceArea(surface_areas[i], "cm2/g");

        ChemicalSystem system(editor);

        ReactionSystem reactions(editor);

        Partition partition(system);
        partition.setKineticSpecies({"Calcite"});

        EquilibriumProblem problem(partition);
        problem.setTemperature(60, "celsius");
        problem.setPressure(100, "bar");
        problem.add("H2O",   1.0, "kg");
        problem.add("O2",    1.0, "umol");
        problem.add("NaCl",  0.9, "mol");
        problem.add("MgCl2", 0.05, "mol");
        problem.add("CaCl2", 0.01, "mol");
        problem.add("CO2",   0.75, "mol");

        // Create initial state by equilibration
        chemical_states.emplace_back(ChemicalState(equilibrate(problem)));
        // Set amount of Calcite
        chemical_states[i].setSpeciesMass("Calcite", 100, "g");

        chemical_states[i].scalePhaseVolume("Aqueous", 1, "m3");
        chemical_states[i].scalePhaseVolume("Calcite", 1, "m3");

        KineticPath path(reactions, partition);

        ChemicalOutput output = path.output();
        output.add("t");
        output.add("pH");
        output.add("speciesMass(Calcite units=g)", "Calcite");
        output.add("speciesMolality(HCO3-)", "HCO3- [molal]");
        output.add("speciesMolality(Ca++)", "Ca++  [molal]");

        std::cout << "**********************************************************************************" << std::endl;
        std::cout << "Surface area: " << surface_areas[i] << " cm2/g" << std::endl;

        output.filename(filename + "-" +  std::to_string(surface_areas[i]) + ".txt");

        double initial_calcite = chemical_states[i].speciesAmount("Calcite");
        path.solve(chemical_states[i], params.t0, params.tfinal);
        double result_calcite = chemical_states[i].speciesAmount("Calcite");
        double dissolved_calcite = initial_calcite - result_calcite;
        std::cout << "Dissol. Calcite: " << dissolved_calcite << " mol" << std::endl;

    }

}


