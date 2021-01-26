#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <ThermoFun/ThermoFun.h>

int main()
{
    ThermoFun::Database database("databases/thermofun/aq17-thermofun.json");

    ChemicalEditor editor(database);
    editor.setTemperatures({500}, "celsius");
    editor.setPressures({3000}, "bar");
    editor.addAqueousPhase({
        "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2",
        "Ca+2", "CaCO3@", "CaCl+", "CaCl2@", "CaHCO3+", "CaHSiO3+", "CaOH+", "CaSiO3@", "K+", "KAlO2@",
        "KCl@", "KOH@", "KCO3-", "KHCO3@", "Mg+2", "MgCO3@", "MgCl+", "MgCl2@", "MgHCO3+", "MgHSiO3+",
        "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaCO3-", "NaCl@", "NaHCO3@", "NaHSiO3@",
        "NaOH@", "HSiO3-", "SiO2@", "CO@", "CO2@", "CO3-2", "HCO3-", "CH4@", "Cl-",
        "HCl@", "H2@", "O2@", "OH-", "H+", "H2O@"});

    editor.addMineralPhase("Albite");
    editor.addMineralPhase("Andalusite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Corundum");
    editor.addMineralPhase("Diopside");
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Enstatite");
    editor.addMineralPhase("Grossular");
    editor.addMineralPhase("Margarite");
    editor.addMineralPhase("Microcline");
    editor.addMineralPhase("Muscovite");
    editor.addMineralPhase("Pargasite-Mg");
    editor.addMineralPhase("Phlogopite");
    editor.addMineralPhase("Quartz");
    editor.addMineralPhase("Sanidine");
    editor.addMineralPhase("Sillimanite");
    editor.addMineralPhase("Zoisite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.add("H2O",             	1000,	"g");
    problem.add("CO2",             	0.001,  "g");
    problem.add("CaCO3",           	1,	    "g");
    problem.add("MgSiO3",          	1,	    "g");
    problem.add("NaCl",            	5,      "g");
    problem.add("NaAlSi3O8",        37,	    "g");
    problem.add("KAl3Si3O10(OH)2",  13,	    "g");
    problem.add("SiO2",          	30,	    "g");
    problem.add("KAlSi3O8",        	20,	    "g");
    problem.setTemperature(500, "celsius");
    problem.setPressure(3000, "bar");

    ChemicalState state = equilibrate(problem);

    state.output("result.txt");
}
