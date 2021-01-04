#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

void MassImbalanceTest()
{
    try
    {
        Database db("supcrt98.xml");
        ChemicalEditor editor(db);
        editor.addAqueousPhaseWithElements("H O Na Cl");
        editor.addMineralPhase("Halite");

        auto system = Reaktoro::ChemicalSystem(editor);
        std::cout << "system = " << system << std::endl;

        auto solver = new Reaktoro::SmartEquilibriumSolver(system);


        Reaktoro::SmartEquilibriumOptions smartOptions;
        smartOptions.reltol = 1;
        smartOptions.learning.nonlinear.max_iterations *= 2;
        smartOptions.learning.optimum.max_iterations *= 2;
        solver->setOptions(smartOptions);

        auto pHMethod = Reaktoro::ChemicalProperty::pH(system);

        double allinput[2][5] = { { 0.144211037489795635, 0.474387123117082921,  0.144209657503628602 , 0.237193561241422823, -1.37935192986337315e-06 },
                                  { 0.144659843404499999, 0.473734453636464936,  0.14469915961447144 , 0.236867226501872669, 3.93168426909351674e-05 }
        };
        double T[2] = { 300, 300 };
        double p[2] = { 32719962.3359830678 , 2.35648e+07 };

        for (int i = 0; i < 2; i += 1)
        {
            Reaktoro::ChemicalState state(system);
            state.setPressure(p[i]);
            state.setTemperature(T[i]);

            Vector input(5);
            for (int j = 0; j < 5; ++j) input[j] = allinput[i][j];

            state.output("solve" + std::to_string(i) + "_in.txt");
            auto res = solver->solve(state, T[i], p[i], input);

            std::cout << std::endl << "Run: " << i << std::endl;
            std::cout << "time: " << res.timing.solve << std::endl;
            std::cout << "smart: " << res.estimate.accepted << std::endl;
            std::cout << "converged: " << res.learning.gibbs_energy_minimization.optimum.succeeded << std::endl;
            std::cout << "Iter: " << res.learning.gibbs_energy_minimization.optimum.iterations << std::endl;
            std::cout << "error: " << res.learning.gibbs_energy_minimization.optimum.error << std::endl;
            std::cout << "Ph: " << double(pHMethod(state.properties())) << std::endl;

            auto elements = state.elementAmounts();
            for (int j=0; j<5; ++j)
            {
                fprintf(stdout, "%s: in: %g, out: %g, proportion: %g\n", system.element(j).name().c_str(), input[j], elements[j], elements[j]/input[j]);
            }

            //fprintf(stdout, "Amount of C: %g\n", state.elementAmount("C"));
            std::cout << "volume: " << state.properties().volume() << std::endl;
            state.output("solve" + std::to_string(i) + ".txt", 18);
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "exception in solve test: " << e.what() << std::endl;
    }
}

void MassImbalanceTest2()
{
    try
    {
        Database db("supcrt98.xml");
        ChemicalEditor editor(db);
        editor.addAqueousPhaseWithElements("H O Na Cl");
        editor.addMineralPhase("Halite");

        auto system = Reaktoro::ChemicalSystem(editor);

        auto solver = new Reaktoro::SmartEquilibriumSolver(system);


        Reaktoro::SmartEquilibriumOptions smartOptions;
        smartOptions.reltol = 1;
        smartOptions.learning.nonlinear.max_iterations *= 2;
        smartOptions.learning.optimum.max_iterations *= 2;
        solver->setOptions(smartOptions);

        auto pHMethod = Reaktoro::ChemicalProperty::pH(system);

        double allinput[2][5] = { { 0.14417546353022187, 0.474425201430425547 ,  0.144181098768655536  , 0.237212600398161805  , 5.6358725353359643e-06 },
                                  { 0.144121875986094744, 0.47451282206401979,  0.144115383293366012 , 0.237256410714771382, -6.49205825196948292e-06 }
        };
        double T[2] = { 300, 300 };
        double p[2] = { 34468612.1448629647 , 34988313.7764134407 };

        for (int i = 0; i < 2; i += 1)
        {
            Reaktoro::ChemicalState state(system);
            state.setPressure(p[i]);
            state.setTemperature(T[i]);

            Vector input(5);
            for (int j = 0; j < 5; ++j) input[j] = allinput[i][j];

            state.output("solve" + std::to_string(i) + "_in.txt");
            auto res = solver->solve(state, T[i], p[i], input);

            std::cout << std::endl << "Run: " << i << std::endl;
            std::cout << "time: " << res.timing.solve << std::endl;
            std::cout << "smart: " << res.estimate.accepted << std::endl;
            std::cout << "converged: " << res.learning.gibbs_energy_minimization.optimum.succeeded << std::endl;
            std::cout << "Iter: " << res.learning.gibbs_energy_minimization.optimum.iterations << std::endl;
            std::cout << "error: " << res.learning.gibbs_energy_minimization.optimum.error << std::endl;
            std::cout << "Ph: " << double(pHMethod(state.properties())) << std::endl;

            auto elements = state.elementAmounts();
            for (int j = 0; j < 5; ++j)
            {
                fprintf(stdout, "%s: in: %.18g, out: %.18g, proportion: %.10g\n", system.element(j).name().c_str(), input[j], elements[j], elements[j] / input[j]);
            }

            //fprintf(stdout, "Amount of C: %g\n", state.elementAmount("C"));
            std::cout << "volume: " << state.properties().volume() << std::endl;
            state.output("solve" + std::to_string(i) + ".txt", 18);
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "exception in solve test: " << e.what() << std::endl;
    }
}
int main()
{
    MassImbalanceTest();
    //MassImbalanceTest2();
}
