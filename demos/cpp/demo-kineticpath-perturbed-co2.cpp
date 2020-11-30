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

    auto folder = params.makeResultsFolder();

    ChemicalEditor editor;
    editor.addAqueousPhaseWithElementsOf("H2O HCl CaCO3");
    editor.addMineralPhase("Calcite");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");
    reaction.setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);

    // Function to evaluate ph
    auto ph_func = ChemicalProperty::pH(system);

    ReactionSystem reactions(editor);

    Partition partition(system);
    partition.setKineticSpecies({"Calcite"});

    EquilibriumProblem problem(partition);
    problem.setTemperature(60, "celsius");
    problem.setPressure(100, "bar");
    problem.add("H2O", 1, "kg");

    // Create a vector of chemical states
    std::vector<ChemicalState> chemical_states;
    // Initialize amount of equilibrium simulation with different C02(g) amounts
    const unsigned int states_num = 5;

    for(unsigned i = 0; i < states_num; ++i){
        // Add CO2 (increases the initial amount by 0.1 mol)
        problem.add("CO2", 0.1, "mol");
        // Create initial state by equilibration
        chemical_states.emplace_back(ChemicalState(equilibrate(problem)));
        // Set amount of Calcite
        chemical_states[i].setSpeciesMass("Calcite", 100, "g");
    }

    KineticPath path(reactions, partition);

    ChemicalOutput output = path.output();
    output.add("t");
    output.add("pH");
    output.add("speciesMass(Calcite units=g)", "Calcite");
    output.add("speciesMolality(HCO3-)", "HCO3- [molal]");
    output.add("speciesMolality(Ca++)", "Ca++  [molal]");

    for(unsigned i = 0; i < states_num; ++i){

        std::cout << "**********************************************************************************" << std::endl;
        std::cout << i + 1 << std::endl;
        std::cout << "**********************************************************************************" << std::endl;

        output.filename(folder + "/path" + std::to_string(i + 1));

        // Output initial values
        double initial_calcite = chemical_states[i].speciesAmount("Calcite");
        double initial_co2 = chemical_states[i].speciesAmount("CO2(aq)");
        double initial_ph = ph_func(chemical_states[i].properties()).val;

        std::cout << "initial Calcite : " << initial_calcite << std::endl;
        std::cout << "initial CO2     : " << initial_co2 << std::endl;
        std::cout << "initial pH      : " << initial_ph << std::endl;

        // Perform kinetic path simulations
        path.solve(chemical_states[i], params.t0, params.tfinal);

        // Output resulting values
        double result_ph = ph_func(chemical_states[i].properties()).val;
        double result_calcite = chemical_states[i].speciesAmount("Calcite");
        double dissolved_calcite = initial_calcite - result_calcite;

        std::cout << "----------------------------" << std::endl;
        std::cout << "result  pH      : " << result_ph << std::endl;
        std::cout << "result  Calcite : " << result_calcite << std::endl;
        std::cout << "dissol. Calcite : " << dissolved_calcite << std::endl;

    }
}