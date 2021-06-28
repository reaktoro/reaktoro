#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database db("supcrt98.yaml");

    AqueousPhase aqueousphase("H2O(aq) H+ OH- HCO3- CO3-2 CO2(aq)");
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    ChemicalSystem system(phases);
    ChemicalState state(system);
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    state.setSpeciesMass("H2O(aq)", 1.0, "kg");
    state.setSpeciesAmount("CO2(g)", 10.0, "mol");

    EquilibriumSolver solver(system);
    solver.solve(state);

    std::cout << "Finished calculation!" << std::endl;
}
