// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
//   • Allan Leal (21 September 2022)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Note: ODML is currently implemented with only an on-demand clustering
// strategy for storing and searching learning data. However, this is not
// suitable for the example below, which would possibly benefit more from
// a nearest-neighbor search algorithm instead, to be implemented.
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <reaktplot/reaktplot.hpp>
using namespace reaktplot;

int main(int argc, char const *argv[])
{
    // -----------------------------------------------------------------------------
    // CONFIGURING CHEMICAL SYSTEM FOR THE CALCULATIONS
    // -----------------------------------------------------------------------------

    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution("H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-2 CO2(aq)");
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gases("CO2(g) H2O(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    ChemicalSystem system(db, solution, gases);

    // -----------------------------------------------------------------------------
    // PARAMETERS FOR THE CALCULATIONS
    // -----------------------------------------------------------------------------

    auto N = 50;  // number of temperature and pressure points

    ArrayXd temperatures = linspace(25.0, 90.0, N);  // temperatures from 25 to 90 °C
    ArrayXd pressures = linspace(1.0, 300.0, N);  // pressures from 1 to 300 bar

    // -----------------------------------------------------------------------------
    // CONFIGURATION OF THE SMART AND EXACT EQUILIBRIUM SOLVERS
    // -----------------------------------------------------------------------------

    SmartEquilibriumSolver smart_solver(system);
    EquilibriumSolver exact_solver(system);

    EquilibriumOptions exact_solver_options;

    SmartEquilibriumOptions smart_solver_options;
    smart_solver_options.reltol = 0.001;

    exact_solver.setOptions(exact_solver_options);
    smart_solver.setOptions(smart_solver_options);

    // The initial chemical state in disequilibrium
    ChemicalState state0(system);
    state0.set("H2O(aq)", 1.0, "kg");
    state0.set("Na+",     1.0, "mol");
    state0.set("Cl-",     1.0, "mol");
    state0.set("CO2(g)", 10.0, "mol");

    ArrayXXd exact_mCO2 = zeros(N, N);
    ArrayXXd smart_mCO2 = zeros(N, N);

    Vec<double> xpredicted, ypredicted; // the points (x, y) = (T, P) where prediction occurred
    Vec<double> xlearning, ylearning;   // the points (x, y) = (T, P) where learning occurred

    // -----------------------------------------------------------------------------
    // CALCULATIONS USING SMART SOLVER RELYING ON ON-DEMAND MACHINE LEARNING (ODML)
    // -----------------------------------------------------------------------------

    std::cout << "Solving equilibrium problems using EquilibriumSolverODML...\n";

    Stopwatch stopwatch;

    ChemicalState state(state0);

    for(auto [i, T] : enumerate(temperatures))
    {
        for(auto [j, P] : enumerate(pressures))
        {
            state.temperature(T, "celsius");
            state.pressure(P, "bar");

            stopwatch.start();

            auto result = smart_solver.solve(state);

            stopwatch.pause();

            errorif(result.failed(), "Equilibrium calculation (with ODML) failed at ", T, " celsius and ", P, " bar.");

            smart_mCO2(j, i) = state.props().elementAmountInPhase("C", "AqueousPhase");

            if(result.predicted())
            {
                xpredicted.push_back(T);
                ypredicted.push_back(P);
            }
            else
            {
                xlearning.push_back(T);
                ylearning.push_back(P);
            }
        }
    }

    auto computing_time_smart_solver = stopwatch.time();

    std::cout << "Finished! Computing time using SmartEquilibriumSolver: " << computing_time_smart_solver << " s\n";

    // -----------------------------------------------------------------------------
    // CALCULATIONS USING CONVENTIONAL/EXACT SOLVER
    // -----------------------------------------------------------------------------

    std::cout << "Solving equilibrium problems using EquilibriumSolver...\n";

    stopwatch.reset();

    state = state0;

    for(auto [i, T] : enumerate(temperatures))
    {
        for(auto [j, P] : enumerate(pressures))
        {
            // print(f"...Progress: {(i*N + j)/(N*N) * 100.0} %")

            state.temperature(T, "celsius");
            state.pressure(P, "bar");

            stopwatch.start();

            auto result = exact_solver.solve(state);

            stopwatch.pause();

            errorif(result.failed(), "Equilibrium calculation (without ODML) failed at ", T, " celsius and ", P, " bar.");

            exact_mCO2(j, i) = state.props().elementAmountInPhase("C", "AqueousPhase");
        }
    }

    auto computing_time_exact_solver = stopwatch.time();

    std::cout << "Finished! Computing time using EquilibriumSolver: " << computing_time_exact_solver << " s\n";

    // -----------------------------------------------------------------------------
    // PLOTTING THE RESULTS
    // -----------------------------------------------------------------------------

    auto speedup = computing_time_exact_solver / computing_time_smart_solver;
    auto prediction_success_rate = double(xpredicted.size())/(N*N) * 100.0;  // in percent

    Figure fig1;
    fig1.title("CO2 SOLUBILITY IN BRINE - WITH ODML - SPEEDUP: " + str(speedup) + " PREDICTED: " + str(prediction_success_rate) + "%");
    fig1.xaxisTitle("Temperature [°C]");
    fig1.yaxisTitle("Pressure [bar]");
    fig1.drawContour(temperatures, pressures, smart_mCO2, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)));
    fig1.drawMarkers(xpredicted, ypredicted, "Prediction", MarkerSpecs().color("black").size(5));
    fig1.show();
    fig1.save("ex-smart-equilibrium-co2-solubility-brine-odml-yes.pdf");

    Figure fig2;
    fig2.title("CO2 SOLUBILITY IN BRINE - WITHOUT ODML");
    fig2.xaxisTitle("Temperature [°C]");
    fig2.yaxisTitle("Pressure [bar]");
    fig2.drawContour(temperatures, pressures, exact_mCO2, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)));
    fig2.show();
    fig2.save("ex-smart-equilibrium-co2-solubility-brine-odml-no.pdf");

    ArrayXXd error = abs(exact_mCO2 - smart_mCO2)/abs(exact_mCO2) * 100.0; // error in percent

    Figure fig3;
    fig3.title("CO2 SOLUBILITY IN BRINE - ERROR USING ODML");
    fig3.xaxisTitle("Temperature [°C]");
    fig3.yaxisTitle("Pressure [bar]");
    fig3.drawContour(temperatures, pressures, error, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)));
    fig3.show();
    fig3.save("ex-smart-equilibrium-co2-solubility-brine-odml-error.pdf");

    return 0;
}
