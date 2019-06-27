#include <Reaktoro/Reaktoro.hpp>
#include <fstream>


#include "W:\release\Projects\Reaktoro\demos\cpp\phaseid.hpp"
#include <Reaktoro/Oilphase/oilphase.cpp>

using namespace Reaktoro;

struct Point {
	double temperature; //[ºC]
	double pressure;    //[bar]
};

struct Interval {
	Point initial_point;
	Point final_point;
};

auto def_points(Interval P, int num_div) -> std::vector<Point> {
	double dT{ (P.final_point.temperature - P.initial_point.temperature) / (double)num_div };
	double dP{ (P.final_point.pressure - P.initial_point.pressure) / (double)num_div };
	std::vector<Point> points;
	points.reserve(num_div * num_div + 2 * num_div + 1);
	for (auto i = 0; i <= num_div; i++) {
		
		for (auto j = 0; j <= num_div; j++) {
			//if (i % 2 == 0) {
				points.push_back({ P.initial_point.temperature + (dT * j), P.initial_point.pressure + (dP * i) });
			//}
			//else {
				//points.push_back({ P.final_point.temperature - (dT * j), P.initial_point.pressure + (dP * i) });
			//}

		}
		
		/*
		for (auto j = 0; j <= num_div; j++) {
			points.push_back({ P.initial_point.temperature + (dT*i), P.initial_point.pressure + (dP*j) });
		}
		*/
	}
	return points;
}

auto count_num_of_phase(const ChemicalState& state) -> int {
	auto num_phase{ 0 };
	for (auto phase : state.system().phases()) {
		if (state.phaseAmount(phase.name()) > 1e-7) {
			num_phase++;
		}
	}
	return num_phase;
}

auto output(const ChemicalState& state, const EquilibriumResult& res, std::ofstream& out) -> void {

	const auto& system = state.system();
	auto num_phases = count_num_of_phase(state);

	out << std::setw(30) << state.temperature();
	out << std::setw(30) << state.pressure();
	out << std::setw(30) << res.optimum.succeeded;

	out << std::setw(30) << num_phases;

	for (auto phase : state.system().phases()) {
		out << std::setw(30) << state.phaseAmount(phase.name());
	}

	auto densities = state.properties().phaseMasses().val / state.properties().phaseVolumes().val;

	for (auto i = 0; i < densities.size(); i++) {
		out << std::setw(30) << densities[i];
	}

	for (auto i = 0; i < state.speciesAmounts().size(); i++) {
		out << std::setw(15) << state.speciesAmounts()[i];
	}

	out << std::endl;
}


int main()
{


	//Database db("supcrt98.xml");
	auto db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml");

	// NOTE: Copy & Pasted from ChemicalEditor {
	// The default temperatures for the interpolation of the thermodynamic properties (in units of celsius)
	std::vector<double> temperatures{ 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300 };

	// The default pressures for the interpolation of the thermodynamic properties (in units of bar)
	std::vector<double> pressures{ 1, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 };

	// Convert the temperatures and pressures to units of kelvin and pascal respectively
	for (auto& x : temperatures) { x = x + 273.15; }
	for (auto& x : pressures) { x = x * 1.0e+5; }

	auto editor = ChemicalEditor(db);

	editor.addAqueousPhase({ "H2O(l)" });

	//// Approach 1:
	//// editor.addGaseousPhase({"H2O(g)", "CH4(g)" /* , "C2H6(g) " */ });
	//editor.addGaseousPhase({"CH4(g)" /* , "C2H6(g) " */ });
	//ChemicalSystem system(editor);


	// Approach 2:
	std::vector<Phase> phases;
	phases.push_back(convertPhase(editor.aqueousPhase(), db));

	{
		auto gas_species = std::vector<GaseousSpecies>{
			//db.gaseousSpecies("H2O(g)"),
			db.gaseousSpecies("CH4(g)"),
			db.gaseousSpecies("CO2(g)"),
			db.gaseousSpecies("H2S(g)"),
			//db.gaseousSpecies("C1(g)"),
			//db.gaseousSpecies("C2(g)"),
			//db.gaseousSpecies("C3(g)"),
			//db.gaseousSpecies("C4(g)"),
			//db.gaseousSpecies("C5(g)"),
			//db.gaseousSpecies("C6(g)"),
			//db.gaseousSpecies("C7(g)"),
		};
		auto mixture = GaseousMixture(gas_species);
		auto gas = GaseousPhase(mixture);

		gas.setChemicalModelSoaveRedlichKwong();
		//gas.setChemicalModelPengRobinson();

		phases.push_back(convertPhase(gas, db));

	}

	{
		auto oil_species = std::vector<HydrocarbonSpecies>{
			HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("H2O(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C1(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C2(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C3(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C4(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C5(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C6(oil)")),
			//HydrocarbonSpecies(db.gaseousSpecies("C7(oil)")),

		};

		//OilSpecies C2H6;
		//C2H6.setName("C2H6(o)");
		//C2H6.setFormula("C2H6");

		//auto C2H6_species = oil_species[0].elements(); // copy from CH4
		//// change values
		//for (auto& kv : C2H6_species) {
		//    if (kv.first.name() == "C") kv.second = 2.0;
		//    if (kv.first.name() == "H") kv.second = 6.0;
		//}

		//// http://www.coolprop.org/fluid_properties/fluids/Ethane.html
		//C2H6.setAcentricFactor(0.099);
		//C2H6.setCriticalPressure(4872200.0); // [Pa]
		//C2H6.setCriticalTemperature(305.322); // [K]
		//C2H6.setElements(C2H6_species);
		//C2H6.setThermoData(oil_species[0].thermoData());

		//oil_species.push_back(C2H6);

		auto mixture = HydrocarbonMixture(oil_species);
		auto oil = HydrocarbonPhase(mixture);

		oil.setChemicalModelSoaveRedlichKwong();

		//oil.setChemicalModelRedlichKwong();
		
		//oil.setChemicalModelPengRobinson();

		phases.push_back(convertPhase(oil, db));
	}


	// auto system = ChemicalSystem(editor);
	auto system = ChemicalSystem(phases);

	auto p = def_points({ {-48.15, 1.0} , {26.85, 130.0 } }, 100);
	//auto p = std::vector<Point>{ {-50.0, 100} };


	/*
	for (auto i = 0; i < 6333; i++) {
		p.erase(p.begin());
	}
	*/
	std::ofstream out("W:\\release\\Projects\\reaktoronote\\notebook\\\oilphase\\pvtlib_validation\\CH4_CO2_H2S\\0_7CH4_0_10H2S_0_20CO2_no_initial_guess.txt");

	out << std::setw(30) << "T" << std::setw(30) << "P" << std::setw(30) << "converged" << std::setw(30) << "num_phases";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_amount";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_density";

	for (auto i = 0; i < system.species().size(); i++) {
		out << std::setw(15) << system.species(i).name();
	}
	out << std::endl;

	//for (auto k = 0; k < 10; k++) {

		auto state = ChemicalState(system);
		/*
		state.setSpeciesAmount("CH4(oil)", 1.0 - (1.0 / 10.0)*k);
		state.setSpeciesAmount("CH4(g)", (1.0 / 10.0)*k);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "CH4(oil)" << 1.0 - (1.0 / 10.0)*k << std::endl;
		std::cout << "CH4(g)" << (1.0 / 10.0)*k << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		*/
		for (auto i = 0; i < p.size(); i++) {

			auto system = ChemicalSystem(phases);

			auto problem = EquilibriumProblem(system);

			problem.setTemperature(p[i].temperature, "celsius");
			problem.setPressure(p[i].pressure, "bar");
			//problem.setTemperature(-100, "celsius");
			//problem.setPressure(10, "bar");
			problem.add("CH4(g)", 0.7, "mol");
			problem.add("H2S(g)", 0.2, "mol");
			problem.add("CO2(g)", 0.1, "mol");
			//problem.add("H2S(g)", 0.888, "mol");
			
			/*
			problem.add("C1(g)", 0.860, "mol");
			problem.add("C2(g)", 0.050, "mol");
			problem.add("C3(g)", 0.050, "mol");
			problem.add("C4(g)", 0.020, "mol");
			problem.add("C5(g)", 0.010, "mol");
			problem.add("C6(g)", 0.005, "mol");
			problem.add("C7(g)", 0.0005, "mol");
			problem.add("C1(oil)", 0.450, "mol");
			problem.add("C2(oil)", 0.050, "mol");
			problem.add("C3(oil)", 0.050, "mol");
			problem.add("C4(oil)", 0.030, "mol");
			problem.add("C5(oil)", 0.010, "mol");
			problem.add("C6(oil)", 0.010, "mol");
			problem.add("C7(oil)", 0.400, "mol");

			state.setSpeciesAmount("C1(g)", 0.860);
			state.setSpeciesAmount("C2(g)", 0.050);
			state.setSpeciesAmount("C3(g)", 0.050);
			state.setSpeciesAmount("C4(g)", 0.020);
			state.setSpeciesAmount("C5(g)", 0.010);
			state.setSpeciesAmount("C6(g)", 0.005);
			state.setSpeciesAmount("C7(g)", 0.005);
			state.setSpeciesAmount("C1(oil)", 0.450);
			state.setSpeciesAmount("C2(oil)", 0.050);
			state.setSpeciesAmount("C3(oil)", 0.050);
			state.setSpeciesAmount("C4(oil)", 0.030);
			state.setSpeciesAmount("C5(oil)", 0.010);
			state.setSpeciesAmount("C6(oil)", 0.010);
			state.setSpeciesAmount("C7(oil)", 0.400);
			*/

			auto solver = EquilibriumSolver(problem.system());
			solver.setPartition(problem.partition());

			auto options = EquilibriumOptions();
			options.hessian = GibbsHessian::Exact;
			options.nonlinear.max_iterations = 100;
			options.optimum.max_iterations = 200;
			options.optimum.output.active = false;
			options.optimum.ipnewton.step = Conservative;
			//options.optimum.tolerance = 1.e-7;
			solver.setOptions(options);
			
			//if (p[i].pressure == 1.0) {
			state.setSpeciesAmounts(0.0);
			
			//(0.925382, 0.074618)
			//(0.118133, 0.881867)
			
			//state.setSpeciesAmount("CH4(g)", 1.0);
			//state.setSpeciesAmount("CO2(g)", 1.0);

			//state.setSpeciesAmount("H2S(g)", 0.062826 );
			
			

			//}


			auto phaseidMethod = phaseIdentificationMethod::Gibbs_residual_based;



			auto res_final = solver.solve(state, problem, phaseidMethod, 0);

			//std::cout << state << std::endl;

			//auto res = solver.solve(state, problem, phaseidMethod, 0);

			//std::cout << state << std::endl;

			output(state, res_final, out);

			//std::cout << state << std::endl;
			//std::cout << res.optimum.succeeded << std::endl;
		}
	//}

	

}