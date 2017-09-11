// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "TestAlgorithmIpOpt.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

auto test_ipfeasible() -> void
{
    const unsigned n = 5;
    const unsigned m = 2;

    OptimumProblem problem(n, m);
    OptimumState state;
    OptimumOptions options;

    ipfeasible(problem, state, options);

    ASSERT_EQUAL(n, state.x.size());
    ASSERT_EQUAL(m, state.y.size());
    ASSERT_EQUAL(n, state.z.size());
    ASSERT_EQUAL(n, state.zu.size());
    ASSERT(arma::all(state.x > 0.0));
    ASSERT(arma::all(state.z == 0.0));
    ASSERT(arma::all(state.zu == 0.0));
}

auto test_ipopt_parabolic() -> void
{
    const unsigned m = 1;
    const unsigned n = 2;

    ObjectiveFunction objective = [](const Vector& x)
    {
        ObjectiveResult f;
        f.func    = x[0]*x[0] + x[1]*x[1];
        f.grad    = 2*x;
        f.hessian = 2*arma::eye(n, n);
        return f;
    };

    ConstraintFunction constraint = [](const Vector& x)
    {
        Matrix A = arma::ones(m, n);
        Vector b = {1.0};
        ConstraintResult h;
        h.func = A*x-b;
        h.grad = A;
        h.hessian = arma::zeros(m, n, n);
        return h;
    };

    OptimumProblem problem(n, m);
    problem.setObjective(objective);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0.0);

    OptimumState state;
    OptimumOptions options;

    ipfeasible(problem, state, options);
    state.x = {2.0, 0.01};
    ipopt(problem, state, options);

    ASSERT_EQUAL_DELTA(0.5, state.x[0], 1e-8);
    ASSERT_EQUAL_DELTA(0.5, state.x[1], 1e-8);
}

auto test_ipopt_logarithmic() -> void
{
    const unsigned m = 1;
    const unsigned n = 2;

    ObjectiveFunction objective = [](const Vector& x)
    {
        const double x0 = x[0];
        const double x1 = x[1];
        ObjectiveResult f;
        f.func    = x0 + x1 + x0*std::log(x0) + x1*std::log(x1);
        f.grad    = 2 + arma::log(x);
        f.hessian = arma::diagmat(1.0/x);
        return f;
    };

    ConstraintFunction constraint = [](const Vector& x)
    {
        Matrix A = arma::ones(m, n);
        Vector b = {1.0};
        ConstraintResult h;
        h.func = A*x-b;
        h.grad = A;
        h.hessian = arma::zeros(m, n, n);
        return h;
    };

    OptimumProblem problem(n, m);
    problem.setObjective(objective);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0.0);

    OptimumState state;
    OptimumOptions options;

    ipfeasible(problem, state, options);
    state.x = {2.0, 0.01};
    ipopt(problem, state, options);

    ASSERT_EQUAL_DELTA(0.5, state.x[0], 1e-15);
    ASSERT_EQUAL_DELTA(0.5, state.x[1], 1e-15);
}

auto createChemicalSystem() -> ChemicalSystem
{
    std::vector<SpeciesThermoModel> aqueous_thermo_models(8);
    aqueous_thermo_models[0].gibbs_energy = [](double, double) { return ThermoScalar(-237.141, 0.0, 0.0); };
    aqueous_thermo_models[1].gibbs_energy = [](double, double) { return ThermoScalar(   0.000, 0.0, 0.0); };
    aqueous_thermo_models[2].gibbs_energy = [](double, double) { return ThermoScalar(-157.297, 0.0, 0.0); };
    aqueous_thermo_models[3].gibbs_energy = [](double, double) { return ThermoScalar(  17.723, 0.0, 0.0); };
    aqueous_thermo_models[4].gibbs_energy = [](double, double) { return ThermoScalar(  16.544, 0.0, 0.0); };
    aqueous_thermo_models[5].gibbs_energy = [](double, double) { return ThermoScalar(-586.940, 0.0, 0.0); };
    aqueous_thermo_models[6].gibbs_energy = [](double, double) { return ThermoScalar(-527.983, 0.0, 0.0); };
    aqueous_thermo_models[7].gibbs_energy = [](double, double) { return ThermoScalar(-385.974, 0.0, 0.0); };

    std::vector<SpeciesThermoModel> gaseous_thermo_models(5);
    gaseous_thermo_models[0].gibbs_energy = [](double, double) { return ThermoScalar(-228.570, 0.0, 0.0); };
    gaseous_thermo_models[1].gibbs_energy = [](double, double) { return ThermoScalar(-394.400, 0.0, 0.0); };
    gaseous_thermo_models[2].gibbs_energy = [](double, double) { return ThermoScalar( -50.659, 0.0, 0.0); };
    gaseous_thermo_models[3].gibbs_energy = [](double, double) { return ThermoScalar(   0.000, 0.0, 0.0); };
    gaseous_thermo_models[4].gibbs_energy = [](double, double) { return ThermoScalar(   0.000, 0.0, 0.0); };

    std::vector<Species> aqueous_species(8);

    aqueous_species[0].setName("H2O(l)");
    aqueous_species[1].setName("H+");
    aqueous_species[2].setName("OH-");
    aqueous_species[3].setName("H2(aq)");
    aqueous_species[4].setName("O2(aq)");
    aqueous_species[5].setName("HCO3-");
    aqueous_species[6].setName("CO3--");
    aqueous_species[7].setName("CO2(aq)");

    aqueous_species[0].setCharge(0);
    aqueous_species[1].setCharge(1);
    aqueous_species[2].setCharge(-1);
    aqueous_species[3].setCharge(0);
    aqueous_species[4].setCharge(0);
    aqueous_species[5].setCharge(-1);
    aqueous_species[6].setCharge(-2);
    aqueous_species[7].setCharge(0);

    aqueous_species[0].setElements({"H", "O"});
    aqueous_species[1].setElements({"H"});
    aqueous_species[2].setElements({"O", "H"});
    aqueous_species[3].setElements({"H"});
    aqueous_species[4].setElements({"O"});
    aqueous_species[5].setElements({"H", "C", "O"});
    aqueous_species[6].setElements({"C", "O"});
    aqueous_species[7].setElements({"C", "O"});

    aqueous_species[0].setElementAtoms({2, 1});
    aqueous_species[1].setElementAtoms({1});
    aqueous_species[2].setElementAtoms({1, 1});
    aqueous_species[3].setElementAtoms({2});
    aqueous_species[4].setElementAtoms({2});
    aqueous_species[5].setElementAtoms({1, 1, 3});
    aqueous_species[6].setElementAtoms({1, 3});
    aqueous_species[7].setElementAtoms({1, 2});

    aqueous_species[0].setThermoModel(aqueous_thermo_models[0]);
    aqueous_species[1].setThermoModel(aqueous_thermo_models[1]);
    aqueous_species[2].setThermoModel(aqueous_thermo_models[2]);
    aqueous_species[3].setThermoModel(aqueous_thermo_models[3]);
    aqueous_species[4].setThermoModel(aqueous_thermo_models[4]);
    aqueous_species[5].setThermoModel(aqueous_thermo_models[5]);
    aqueous_species[6].setThermoModel(aqueous_thermo_models[6]);
    aqueous_species[7].setThermoModel(aqueous_thermo_models[7]);

    std::vector<Species> gaseous_species(5);

    gaseous_species[0].setName("H2O(g)");
    gaseous_species[1].setName("CO2(g)");
    gaseous_species[2].setName("CH4(g)");
    gaseous_species[3].setName("H2(g)");
    gaseous_species[4].setName("O2(g)");

    gaseous_species[0].setCharge(0);
    gaseous_species[1].setCharge(0);
    gaseous_species[2].setCharge(0);
    gaseous_species[3].setCharge(0);
    gaseous_species[4].setCharge(0);

    gaseous_species[0].setElements({"H", "O"});
    gaseous_species[1].setElements({"C", "O"});
    gaseous_species[2].setElements({"C", "H"});
    gaseous_species[3].setElements({"H"});
    gaseous_species[4].setElements({"O"});

    gaseous_species[0].setElementAtoms({2, 1});
    gaseous_species[1].setElementAtoms({1, 2});
    gaseous_species[2].setElementAtoms({1, 4});
    gaseous_species[3].setElementAtoms({2});
    gaseous_species[4].setElementAtoms({2});

    gaseous_species[0].setThermoModel(gaseous_thermo_models[0]);
    gaseous_species[1].setThermoModel(gaseous_thermo_models[1]);
    gaseous_species[2].setThermoModel(gaseous_thermo_models[2]);
    gaseous_species[3].setThermoModel(gaseous_thermo_models[3]);
    gaseous_species[4].setThermoModel(gaseous_thermo_models[4]);

    std::vector<PhaseThermoModel> phase_thermo_models(2);

    phase_thermo_models[0].activity = [](double T, double P, const Vector& n)
    {
        return moleFractions(n);
    };

    phase_thermo_models[1].activity = [](double T, double P, const Vector& n)
    {
        const double Pb = convert<Pa,bar>(P);
        ChemicalVector x = moleFractions(n);
        ChemicalVector a(x.val*Pb, x.ddT()*Pb, x.ddP()*Pb, x.ddn*Pb);
        return a;
    };

    std::vector<Species> aqueous_species2;
    aqueous_species2.push_back(aqueous_species[0]);
    aqueous_species2.push_back(aqueous_species[1]);
    aqueous_species2.push_back(aqueous_species[2]);
    aqueous_species2.push_back(aqueous_species[5]);
    aqueous_species2.push_back(aqueous_species[7]);
    std::vector<Species> gaseous_species2;
    gaseous_species2.push_back(gaseous_species[0]);
    gaseous_species2.push_back(gaseous_species[1]);

    std::vector<Phase> phases(2);
    phases[0].setName("Aqueous");
    phases[0].setSpecies(aqueous_species);
    phases[0].setThermoModel(phase_thermo_models[0]);
    phases[1].setName("Gaseous");
    phases[1].setSpecies(gaseous_species);
    phases[1].setThermoModel(phase_thermo_models[1]);

    ChemicalSystem system(phases);

    return system;
}

auto test_ipopt_equilibrium() -> void
{
    /// The pressure of the system (in units of Pa)
    const double P = 1.0e5;

    /// The temperature of the system (in units of K)
    const double T = 298.15;

    /// The universal gas constant (in units of kJ/(mol*K))
    const double R = 8.3144621e-3;

    const double nH2O = 55.1;
    const double nCO2 = 0.5;

    ChemicalSystem system = createChemicalSystem();

    Matrix A = formulaMatrix(system);

    A.resize(A.n_rows + 1, A.n_cols);
    A.row(A.n_rows-1) = speciesCharges(system.species()).t();

    auto elements = system.elements();

    Vector b = {nCO2, 2*nH2O, nH2O + 2*nCO2, 0.0};

    Vector u0 = standardGibbsEnergies(system, T, P).val;

    ObjectiveFunction objective = [=](const Vector& n)
    {
        ChemicalVector a = activities(system, T, P, n);
        Vector u = u0/(R*T) + arma::log(a.val);
        Matrix dudn = arma::diagmat(1/a.val) * a.ddn;

        ObjectiveResult f;
        f.func    = arma::dot(n, u);
        f.grad    = u;
        f.hessian = dudn;
        return f;
    };

    ConstraintFunction constraint = [&](const Vector& n)
    {
        ConstraintResult h;
        h.func = A*n-b;
        h.grad = A;
        return h;
    };

    const unsigned N = numSpecies(system);
    const unsigned E = numElements(system);

//    const Vector lower = 1e-14 * arma::ones(N);
//    const Vector lower = arma::zeros(N);

    OptimumProblem problem(N, E + 1);
    problem.setObjective(objective);
    problem.setConstraint(constraint);
    problem.setLowerBounds(0.0);

    OptimumState state;
    OptimumOptions options;
//    options.ipopt.eta_phi = 1e-8;
    options.output.active = true;
    options.max_iterations = 500;

    ipfeasible(problem, state, options);
    Vector n = 1e-7*arma::ones(N);
    n[indexSpecies(system, "H2O(l)")] = nH2O;
    n[indexSpecies(system, "CO2(g)")] = nCO2;
    state.x  = n;
    state.y  = arma::zeros(E + 1);
    state.z = arma::ones(N);
//    ipopt(problem, state, options);
//    ipopt(problem, state, options); ASSERT(state.statistics.converged);
    options.ipopt.mu = {1e-8, 1e-10, 1e-50};  ipopt(problem, state, options); ASSERT(state.statistics.converged);
//    ipopt(problem, state, options); ASSERT(state.statistics.converged);

    std::cout << "num_iterations: " << state.statistics.num_iterations << std::endl;

//    ASSERT_EQUAL_DELTA(0.5, state.x[0], 1e-15);
//    ASSERT_EQUAL_DELTA(0.5, state.x[1], 1e-15);
}

} // namespace

auto testSuiteAlgorithmIpOpt() -> cute::suite
{
    cute::suite s;

//    s += CUTE(test_ipfeasible);
//    s += CUTE(test_ipopt_parabolic);
//    s += CUTE(test_ipopt_logarithmic);
//    s += CUTE(test_ipopt_equilibrium);

    return s;
}

} // namespace Reaktoro
