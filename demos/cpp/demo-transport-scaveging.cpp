// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <Reaktoro/Reaktoro.hpp>
#include <sstream>
#include <sys/stat.h>

// for windows mkdir
#ifdef _WIN32
#include <direct.h>
#endif

using namespace Reaktoro;

int main()
{
    auto path = "results_demo_transport_and_scaveging"; //folder 
#ifdef _WIN32
    return ::_mkdir(path);
#else
    return ::mkdir(path, 0755);
#endif
 
    auto second = 1;
    auto minute = 60 * second;
    auto hour = 60 * minute;
    auto day = 24 * hour;
    auto year = 365 * day;

    
    //Parameters for the reactive transport simulation
    auto nsteps = 100;      // the number of steps in the reactive transport simulation
    auto ncells = 100;       // the number of cells in the discretization
    auto xl = 0.0;           // the x - coordinate of the left boundary
    auto xr = 100.0;         // the x - coordinate of the right boundary
    auto D = 0.0;            // the diffusion coefficient(in units of m2 / s)
    auto v = 1.05e-5;        // the fluid pore velocity(in units of m / s)
    auto dt = 0.0416*day;    // the time step(in units of s)
    auto T = 25.0;           // the temperature(in units of degC)
    auto P = 1.01325;        // the pressure(in units of bar)

     Database database("supcrt07.xml");
 
     DebyeHuckelParams dhModel{};
     dhModel.setPHREEQC();

     ChemicalEditor editor(database);
     
     editor.addAqueousPhase({ "H2O(l)",  "H+", "OH-", 
                            "HCO3-", "Mg(HCO3)+", "Ca(HCO3)+", "MgCO3(aq)",  "CO3--", "CaCO3(aq)" ,
                            "Ca++", "CaSO4(aq)", "CaOH+", 
                            "Cl-", "FeCl++", "FeCl2(aq)", "FeCl+", 
                            "Fe++", "FeOH+",  "FeOH++", "Fe+++", 
                            "H2(aq)",
                            "K+", "KSO4-", 
                            "Mg++", "MgSO4(aq)", "MgCO3(aq)", "MgOH+", 
                            "Na+", "NaSO4-",
                            "O2(aq)",
                            "H2S(aq)", "HS-", "S5--", "S4--", "S3--", "S2--",
                            "SO4--", "NaSO4-", "MgSO4(aq)", "CaSO4(aq)", "KSO4-", "HSO4-"}).setChemicalModelDebyeHuckel(dhModel);

     editor.addMineralPhase("Pyrrhotite");
     editor.addMineralPhase("Siderite");

     ChemicalSystem system(editor);
 
     EquilibriumInverseProblem problem_ic(system);
     problem_ic.setTemperature(T, "celsius");
     problem_ic.setPressure(P, "bar");
     problem_ic.add("H2O", 58.0, "kg");
     problem_ic.add("Cl-", 1122.3e-3, "kg");
     problem_ic.add("Na+", 624.08e-3, "kg");
     problem_ic.add("SO4--", 157.18e-3, "kg");
     problem_ic.add("Mg++", 74.820e-3, "kg");
     problem_ic.add("Ca++", 23.838e-3, "kg");
     problem_ic.add("K+", 23.142e-3, "kg");
     problem_ic.add("HCO3-", 8.236e-3, "kg");
     problem_ic.add("O2(aq)", 58e-12, "kg");
     problem_ic.add("Pyrrhotite", 0.0, "mol");
     problem_ic.add("Siderite", 0.5, "mol");
     problem_ic.pH(8.951);
     problem_ic.pE(8.676);

 
     EquilibriumInverseProblem problem_bc(system);
     problem_bc.setTemperature(T, "celsius");
     problem_bc.setPressure(P, "bar");
     problem_bc.add("H2O", 58.0, "kg");
     problem_bc.add("Cl-", 1122.3e-3, "kg");
     problem_bc.add("Na+", 624.08e-3, "kg");
     problem_bc.add("SO4--", 157.18e-3, "kg");
     problem_bc.add("Mg++", 74.820e-3, "kg");
     problem_bc.add("Ca++", 23.838e-3, "kg");
     problem_bc.add("K+", 23.142e-3, "kg");
     problem_bc.add("HCO3-", 8.236e-3, "kg");
     problem_bc.add("O2(aq)", 58e-12, "kg");
     problem_bc.add("Pyrrhotite", 0.0, "mol");
     problem_bc.add("Siderite", 0.0, "mol");
     problem_bc.add("HS-", 0.0196504, "mol");
     problem_bc.add("H2S(aq)", 0.167794, "mol");
     problem_bc.pH(5.726);
     problem_bc.pE(8.220);

     ChemicalState state_ic = equilibrate(problem_ic);
     ChemicalState state_bc = equilibrate(problem_bc);

     Mesh mesh(ncells, xl, xr);
 
     ChemicalField field(mesh.numCells(), state_ic);
 
     ReactiveTransportSolver rt(system);
     rt.setMesh(mesh);
     rt.setVelocity(v);
     rt.setDiffusionCoeff(D);
     rt.setBoundaryState(state_bc);
     rt.setTimeStep(dt);
 
 
     auto output = rt.output();
     output.filename(path+"\\reative_transport_siderite_pyrrhotite_pyrite.txt");
     output.add("pH");
     output.add("speciesMolality(H+)");
     output.add("speciesMolality(HS-)");
     output.add("speciesMolality(S2--)");
     output.add("speciesMolality(SO4--)");
     output.add("speciesMolality(HSO4-)");
     output.add("speciesMolality(H2S(aq))");
     output.add("phaseAmount(Pyrrhotite)");
     output.add("phaseAmount(Siderite)");
 
     rt.initialize(field);
 
     auto t = 0.0;
     auto step = 0.0;
 
     while (step <= nsteps) {
         std::cout << step << std::endl;
    
         rt.step(field);
 
         t += dt;
         step += 1;
     }

}