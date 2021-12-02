// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>


// C++ includes
#include<fstream>

// Reaktoro includes
#include <Reaktoro/Transport/ReactiveTransportSolver.hpp>
#include <Reaktoro/Transport/ChemicalField.hpp>
#include <Reaktoro/Transport/ReactiveTransportOptions.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDrummond.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

using namespace Reaktoro;

TEST_CASE("Testing ReactiveTransportSolver", "[ReactiveTransportSolver]")
{
    // Define the list of selected species
    StringList selected_species = "H2O(aq) H+ OH- Na+ Cl- NaCl(aq) Ca+2 Mg+2 HCO3- CO2(aq) CO3-2 CaCl+ CaCl2(aq) O2(aq) Ca(HCO3)+ MgCl+ Mg(HCO3)+";

    // Initialize a thermodynamic database
    SupcrtDatabase db("supcrtbl");

    // Create an aqueous phase automatically selecting all species with given elements, excluding species with tag `organic`
    AqueousPhase aqueousphase(selected_species);

    WHEN("HFK model is used")
    {

        // HKF full system
        aqueousphase.setActivityModel(chain(ActivityModelHKF(), ActivityModelDrummond("CO2")));

        // Create a mineral phase
        MineralPhases mineralphases("Quartz Calcite Dolomite");

        // Construct the chemical system
        ChemicalSystem system(db, aqueousphase, mineralphases);

        // Define equilibrium solver and equilibrate given initial state
        EquilibriumSolver solver(system);

        // Step **: Create the Partition of inert and equilibrium species
        Partition partition(system);

        // Define aqueous and chemical properties
        AqueousProps aprops(system);
        ChemicalProps props(system);

        // Define initial equilibrium state
        ChemicalState state_ic(system);
        state_ic.setTemperature(60.0, "celsius");
        state_ic.setPressure(100.0, "bar");
        state_ic.setSpeciesMass("H2O(aq)", 1.00, "kg");
        state_ic.setSpeciesAmount("O2(aq)", 1.00, "umol");
        state_ic.setSpeciesAmount("NaCl(aq)", 0.70, "mol");
        state_ic.setSpeciesAmount("MgCl+", 1e-10, "mol");
        state_ic.setSpeciesAmount("Cl-", 1e-10, "mol");
        state_ic.setSpeciesAmount("Calcite", 10.00, "mol");
        state_ic.setSpeciesAmount("Quartz", 10.00, "mol");

        solver.solve(state_ic);
        aprops.update(state_ic);

        CHECK(aprops.pH() == Approx(9.26362));

        // Define boundary equilibrium state
        ChemicalState state_bc(system);
        state_bc.setTemperature(60.0, "celsius");
        state_bc.setPressure(100.0, "bar");
        state_bc.setSpeciesMass("H2O(aq)", 1.00, "kg");
        state_bc.setSpeciesAmount("O2(aq)", 1.00, "umol");
        state_bc.setSpeciesAmount("NaCl(aq)", 0.90, "mol");
        state_bc.setSpeciesAmount("MgCl+", 0.05, "mol");
        state_bc.setSpeciesAmount("Cl-", 0.05, "mol");
        state_bc.setSpeciesAmount("CaCl2(aq)", 0.01, "mol");
        state_bc.setSpeciesAmount("CO2(aq)", 0.75, "mol");

        solver.solve(state_bc);
        aprops.update(state_bc);
        CHECK(aprops.pH() == Approx(3.11655));

        // Step **: Scale the boundary condition state
        state_bc.scaleVolume(1.0, "m3");

        // Step **: Scale the volumes of the phases in the initial condition
        state_ic.scalePhaseVolume("AqueousPhase", 0.1, "m3");    // 10% if the 1.0m3
        state_ic.scalePhaseVolume("Quartz", 0.882,
                                  "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% quartz of the rock)
        state_ic.scalePhaseVolume("Calcite", 0.018,
                                  "m3");  // 0.018 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% calcite of the rock)

        // Step **: Create the mesh for the column
        Mesh mesh(100, 0.0, 1.0);

        // Step **: Create a chemical field object with every cell having state given by state_ic
        ChemicalField field(mesh.numCells(), state_ic);

        // Step **: Define the options for the reactive transport solver
        ReactiveTransportOptions reactive_transport_options;

        // Step **: Define the reactive transport modeling
        ReactiveTransportSolver rtsolver(partition);
        rtsolver.setOptions(reactive_transport_options);
        rtsolver.setMesh(mesh);
        rtsolver.setVelocity(1.65344e-06);
        rtsolver.setDiffusionCoeff(1e-09);
        rtsolver.setBoundaryState(state_bc);
        rtsolver.setTimeStep(1800);
        rtsolver.initialize();

        // Step **: Define the quantities that should be output for every cell, every time step
        ChemicalOutput output(rtsolver.output());
        output.add("pH");;
        output.add("speciesAmount(H+)");
        output.add("speciesAmount(Ca+2)");
        output.add("speciesAmount(Mg+2)");
        output.add("speciesAmount(CO3-2)");
        output.add("speciesAmount(CaCl+)");
        output.add("speciesAmount(MgCl+)");
        output.add("speciesAmount(Calcite)");
        output.add("speciesAmount(Dolomite)");
        output.filename("step.txt");

        // Perform one reactive transport time step
        rtsolver.step(field);

        // Load result file and corresponding output string for it
        std::ifstream readfile("step-0.txt");
        std::string output_string;

        if(readfile.is_open())
        {
            SECTION("Checking the string description of chemical output values")
            {
                readfile >> output_string;
                CHECK( output_string == "pH" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(H+)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(Ca+2)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(Mg+2)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(CO3-2)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(CaCl+)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(MgCl+)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(Calcite)" );
                readfile >> output_string;
                CHECK( output_string == "speciesAmount(Dolomite)" );
            }

            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(5.37154)); // pH
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.000635601)); // H+
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(2.27551)); // Ca+2
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.410072)); // Mg+2
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.000279756)); // CO3-2
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.207917)); // CaCl+
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.0431896)); // MgCl+
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(484.542)); // Calcite
            readfile >> output_string;
            CHECK(std::stod(output_string) == Approx(0.587539)); // Dolomite

        }
    }
}
