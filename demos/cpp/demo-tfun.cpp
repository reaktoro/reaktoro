#include "Reaktoro/Reaktoro.hpp"
#include "ThermoFun/ThermoFun.h"

using namespace Reaktoro;

int main()
{
    ThermoFun::Database database("databases/thermofun/aq17.json");
    Time start = time();
    ChemicalEditor editor(database);
    editor.setTemperatures({650}, "celsius");
    editor.setPressures({2000}, "bar");
//    editor.addAqueousPhase("H2O Al3O3 K2O Na2O SiO2 NaCl");
    editor.addAqueousPhase({"Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2", "Cl-", "H+", "H2@", "H2O@", "HCl@", "HSiO3-", "K+", "KAlO2@", "KCl@", "KOH@", "Na+", "NaAl(OH)4@", "NaCl@", "NaHSiO3@", "NaOH@", "O2@", "OH-", "SiO2@"});

//    DebyeHuckelParams params;
//    params.bneutral(0.064);
//    params.aion(3.72);
//    params.bion(0.064);
//    editor.aqueousPhase().setChemicalModelDebyeHuckel(params);

//    editor.addGaseousPhase("CO2");
    editor.addMineralPhase("Albite");
    editor.addMineralPhase("Sanidine");
    editor.addMineralPhase("Andalusite");
    editor.addMineralPhase("Quartz");

    ChemicalSystem system(editor);

    double time_initialization = elapsed(start);

    start = time();

    EquilibriumProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("Al2O3", 183.9, "g");
    problem.add("K2O", 34, "g");
    problem.add("Na2O", 67.1, "g");
    problem.add("SiO2", 715.1, "g");
    problem.add("NaCl", 1, "mol");
    problem.setTemperature(650, "celsius");
    problem.setPressure(2000, "bar");

//    auto props = system.properties(650, 2000);

//    std::cout << props.standardPartialMolarGibbsEnergies().val << std::endl;

//    EquilibriumOptions options;
//    options.optimum.output.active = true;
//    ChemicalState stateaux(system);
//    auto properties = stateaux.properties();
//    std::cout << properties.standardPartialMolarGibbsEnergies().val << std::endl;
    //stateaux.output("stateaux.txt");
    ChemicalState state = equilibrate(problem);

    double time_solve = elapsed(start);

    state.output("result.txt");

    start = time();

    waterDensityWagnerPruss(Temperature(650+273), Pressure(2000e5), StateOfMatter::Liquid);

    double time_watereos = elapsed(start);

    std::cout << "time_initialization = " << time_initialization << std::endl;
    std::cout << "time_solve = " << time_solve << std::endl;
    std::cout << "time_watereos = " << time_watereos << std::endl;
}



