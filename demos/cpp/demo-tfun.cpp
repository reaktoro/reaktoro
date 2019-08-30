#include "Reaktoro/Reaktoro.hpp"
#include "ThermoFun/ThermoFun.h"

using namespace Reaktoro;

int main()
{
    Time start = time();
    Time startT = time();
    ThermoFun::Database database("aq17.json");
    double time_parse_db = elapsed(start);

    start = time();
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

    double time_initialization_editor = elapsed(start);

    start = time();

    ChemicalSystem system(editor);

    double time_initialization_system = elapsed(start);

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

    double time_total = elapsed(startT);

    state.output("result.txt");

    start = time();

    waterDensityWagnerPruss(Temperature(650+273), Pressure(2000e5), StateOfMatter::Liquid);

    ThermoFun::ThermoEngine engine(database);
    ThermoFun::ThermoEngine engine2(database);
    ThermoFun::ThermoBatch batch(database);
    batch.setSolventSymbol("H2O@");

    double time_watereos = elapsed(start);

    start = time();
    double P = 2000e5;
    engine.propertiesSolvent(650+273, P, "H2O@");
    double time_watereos_fun = elapsed(start);

    start = time();
    P = 2000e5;
    engine.electroPropertiesSolvent(650+273, P, "H2O@");
    double time_waterelectro_fun = elapsed(start);

    start = time();
    P = 2000e5;
    engine.propertiesSolvent(650+273, P, "H2O@");
    double time_watereos_fun2 = elapsed(start);

 //   batch.thermoPropertiesSubstance(650+273, P, {"H2O@", "Al(OH)3@"},  {"gibbs_energy","entropy", "volume", "enthalpy"});
    start = time();
    P = 2000e5;
    batch.thermoPropertiesSubstance(650+273, P, {"Albite", "Sanidine", "Andalusite", "Quartz", "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2", "Cl-", "H+", "H2@", "H2O@", "HCl@", "HSiO3-", "K+", "KAlO2@", "KCl@", "KOH@", "Na+", "NaAl(OH)4@", "NaCl@", "NaHSiO3@", "NaOH@", "O2@", "OH-", "SiO2@"}, {"gibbs_energy"});

    engine.propertiesSolvent(650+273, P, "H2O@");
    double time_species_fun = elapsed(start);

    start = time();
    P = 2000e5;
    batch.thermoPropertiesSubstance(650+273, P, {"Albite", "Sanidine", "Andalusite", "Quartz", "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2", "Cl-", "H+", "H2@", "H2O@", "HCl@", "HSiO3-", "K+", "KAlO2@", "KCl@", "KOH@", "Na+", "NaAl(OH)4@", "NaCl@", "NaHSiO3@", "NaOH@", "O2@", "OH-", "SiO2@"}, {"gibbs_energy"});

    engine.propertiesSolvent(650+273, P, "H2O@");
    double time_species_fun2 = elapsed(start);

    start =time();

    engine2.thermoPropertiesSubstance(650+273, P, "Al(OH)2+");
    double time_Alspecies_fun = elapsed(start);

    std::cout << "time_parse_db = " << time_parse_db << std::endl;
    std::cout << "time_initialization_editor = " << time_initialization_editor << std::endl;
    std::cout << "time_initialization_system = " << time_initialization_system << std::endl;
    std::cout << "time_solve = " << time_solve << std::endl;
    std::cout << "time_total_parse_init_solve = " << time_total << std::endl;
    std::cout << "time_watereos_reak = " << time_watereos << std::endl;
    std::cout << "time_watereos_fun = " << time_watereos_fun << std::endl;
//    std::cout << "time_waterelectro_fun = " << time_waterelectro_fun << std::endl;
    std::cout << "time_watereos_fun called 2nd time = " << time_watereos_fun2 << std::endl;
    std::cout << "time_species_fun = " << time_species_fun << std::endl;
    std::cout << "time_species_fun called 2nd time = " << time_species_fun2 << std::endl;
//    std::cout << "time_Alspecies_fun = " << time_Alspecies_fun << std::endl;
}



