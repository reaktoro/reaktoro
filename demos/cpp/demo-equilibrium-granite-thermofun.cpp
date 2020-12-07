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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <ThermoFun/ThermoFun.h>

int main()
{
    // Thermodynamic conditions
    double T = 300.0;
    double P_sat = Reaktoro::waterSaturatedPressureWagnerPruss(Temperature(T + 273.15)).val * 1e-5; // is Psat a water saturation at the T = 300

    ThermoFun::Database database("databases/thermofun/aq17-thermofun.json");

    ChemicalEditor editor(database);
    editor.setTemperatures({T}, "celsius");
    editor.setPressures({P_sat}, "bar");

    // Missing: KAl(OH)4@
    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ "
                                  "K+ KOH@ KCl@ KAlO2@ Al+3 AlOH+2 Al(OH)2+ Al(OH)3@ Al(OH)4-";

    editor.addAqueousPhase(selected_species); // total 25 species
    editor.addMineralPhase("Quartz"); // SiO2
    editor.addMineralPhase("Diaspore"); // AlO(OH)
    editor.addMineralPhase("Gibbsite"); // Al(OH)3
    editor.addMineralPhase("Andalusite"); // Al2SiO5
    editor.addMineralPhase("Kyanite"); // Al2SiO5
    editor.addMineralPhase("Sillimanite"); // Al2SiO5
    editor.addMineralPhase("Muscovite"); // KAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Paragonite"); // NaAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Pyrophyllite"); // Al2Si4O10(OH)2
    editor.addMineralPhase("Kaolinite"); // Al2Si2O5(OH)4
    editor.addMineralPhase("Albite"); // Na(AlSi3)O8
    editor.addMineralPhase("Microcline"); // Microcline is K-feldspar at low T, T < 600

    ChemicalSystem system(editor);
    std::cout << system << std::endl;

    EquilibriumProblem problem(system);
    problem.setTemperature(T, "celsius");
    problem.setPressure(P_sat, "atm");
    // Residual fluid
    problem.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem.add("NaCl", 0.27, "mol"); // NaCl (aq) 0.27 M
    problem.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M

    ChemicalState state = equilibrate(problem);

    state.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state.scalePhaseVolume("Quartz", 0.3 * 0.9, "m3");    // 30% of 90% of remaining volume
    state.scalePhaseVolume("Muscovite", 0.05 * 0.9, "m3");    // 5% of 90% of remaining volume
    state.scalePhaseVolume("Albite", 0.33 * 0.9, "m3");    // 33% of 90% of remaining volume
    state.scalePhaseVolume("Microcline", 0.32 * 0.9, "m3");    // 32% of 90% of remaining volume

    std::cout << state << std::endl;
}
