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
				//points.push_back({ P.initial_point.temperature + (dT * j), P.initial_point.pressure + (dP * i) });
			//}
			//else {
				//points.push_back({ P.final_point.temperature - (dT * j), P.initial_point.pressure + (dP * i) });
			//}

		}
		
		
		for (auto j = 0; j <= num_div; j++) {
			points.push_back({ P.initial_point.temperature + (dT*i), P.final_point.pressure - (dP*j) });
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

auto build_flash(ChemicalSystem& system, std::vector<Point> points, const std::string& file_path) -> void {

	std::ofstream out(file_path);



	out << std::setw(30) << "T" << std::setw(30) << "P" << std::setw(30) << "converged" << std::setw(30) << "num_phases";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_amount";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_density";

	for (auto i = 0; i < system.species().size(); i++) {
		out << std::setw(15) << system.species(i).name();
	}
	out << std::endl;

	for (auto i = 0; i < points.size(); i++) {

		auto problem = EquilibriumProblem(system);

		problem.setTemperature(points[i].temperature, "degC");
		problem.setPressure(points[i].pressure, "bar");
		problem.add("CO2(g)", 0.05, "mol");
		problem.add("H2S(g)", 0.40, "mol");
		problem.add("H2O(aq)", 0.50, "mol");
		//problem.add("O2", 1, "umol");
		//problem.add("H2O", 1, "kg");
		problem.add("CH4(oil)", 0.05, "mol");
		//problem.add("H2S(oil)", 1, "mol");

		//const auto& system = problem.system();

		problem.setTemperature(points[i].temperature, "degC");
		problem.setPressure(points[i].pressure, "bar");

		auto state = ChemicalState(problem.system());
		state.setTemperature(points[i].temperature, "degC");
		state.setPressure(points[i].pressure, "bar");

		auto solver = EquilibriumSolver(problem.system());
		solver.setPartition(problem.partition());

		auto options = EquilibriumOptions();
		options.hessian = GibbsHessian::Exact;
		options.nonlinear.max_iterations = 100;
		options.optimum.max_iterations = 200;

		solver.setOptions(options);

		state.setSpeciesAmounts(0.0);

		auto phaseidMethod = phaseIdentificationMethod::Gibbs_residual_based;
		auto res = solver.solve(state, problem, phaseidMethod, 0);

		/*
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::VolumeMethod;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::workanalisysPengRobinson;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::CriticalPointMethods;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::workanalisysPengRobinson;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod);
		}
		*/
		output(state, res, out);

		std::cout << state << std::endl;

		if (res.optimum.succeeded == false) {
			state.setSpeciesAmounts(0.0);
		}
	}

}

double dist(Point p1, Point p2) {
	return std::sqrt(std::pow(p1.temperature - p2.temperature, 2) + std::pow(p1.pressure - p2.pressure, 2));
}

void distance_order(std::vector<Point>& p) {
	auto min_distance = 1000000000000;
	auto position_of_min_distante = p.size();
	for (auto i = 0; i < p.size(); i++) {
		for (auto j = i + 1; j < p.size(); j++) {
			auto distance = dist(p[i], p[j]);
			if (distance < min_distance) {
				min_distance = distance;
				position_of_min_distante = j;
			}
		}
		auto temp = p[i + 1];
		p[i + 1] = p[position_of_min_distante];
		p[position_of_min_distante] = temp;
		min_distance = 1000000000000;
		position_of_min_distante = p.size();
	}
}

void run_mixture_Huang_1984() {

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
	//editor.addAqueousPhase({ "H2O(l)", "CO2(aq)", "H2S(aq)", "Methane(aq)" });
	editor.addAqueousPhase({ "H2O(l)", "CO2(aq)", "H2S(aq)", "Methane(aq)" });


	// Approach 2:
	std::vector<Phase> phases;
	phases.push_back(convertPhase(editor.aqueousPhase(), db));

	{
		auto gas_species = std::vector<GaseousSpecies>{
			//db.gaseousSpecies("C2H4(g)"),
			db.gaseousSpecies("CO2(g)"),
			db.gaseousSpecies("H2O(g)"),
			db.gaseousSpecies("H2S(g)"),
			//db.gaseousSpecies("O2(g)"),
			db.gaseousSpecies("CH4(g)"),
		};
		auto mixture = GaseousMixture(gas_species);
		auto gas = GaseousPhase(mixture);
		//gas.setChemicalModelRedlichKwong();

		phases.push_back(convertPhase(gas, db));
	}

	{
		auto oil_species = std::vector<HydrocarbonSpecies>{
			HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")),
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
		auto oil = OilPhase(mixture);

		//oil.setChemicalModelRedlichKwong();

		phases.push_back(convertPhase(oil, db));
	}


	// auto system = ChemicalSystem(editor);
	auto system = ChemicalSystem(phases);


	/*
	auto problem = EquilibriumProblem(system);

	problem.setTemperature(79, "degC");
	problem.setPressure(1, "bar");
	problem.add("CO2(g)", 0.05, "mol");
	problem.add("H2S(g)", 0.40, "mol");
	problem.add("H2O(aq)", 0.50, "mol");
	//problem.add("O2", 1, "umol");
	//problem.add("H2O", 1, "kg");
	problem.add("CH4(oil)", 0.05, "mol");
	//problem.add("H2S(oil)", 1, "mol");
	*/
	/*
	{
		auto state = ChemicalState(system);

		auto solver = EquilibriumSolver(system);

		auto res = solver.solve(state, problem);

		std::cout << state << std::endl;

		std::cout << res.optimum.succeeded << std::endl;
	}
	*/

	auto p = def_points({ {0, 52.1},{75,120} }, 10);
	//auto p = def_points({ {0, 0},{10,10} }, 10);

	build_flash(system, p, "W:\\release\\Projects\\reaktoronote\\notebook\\oilphase\\Z_selection\\phase_diagram\\mixture2\\phasediagram_mix2_vol1_phase_id_gibbs_onlyaffter_50_guess_no_initialguess.txt");

}

void phaseBehavior_given_temp(ChemicalSystem system, std::vector<double> pressure_range, const std::string& file_path) {
	std::ofstream out(file_path);

	//const auto& system = problem.system();
	auto state = ChemicalState(system);

	auto solver = EquilibriumSolver(system);

	auto options = EquilibriumOptions();
	options.hessian = GibbsHessian::Exact;
	options.nonlinear.max_iterations = 100;
	options.optimum.max_iterations = 200;

	solver.setOptions(options);

	auto phaseidMethod = phaseIdentificationMethod::workanalisysPengRobinson;

	auto dP = (pressure_range[1] - pressure_range[0]) / 1000;

	auto count = 0;
	auto i = 0;

	std::vector<double> CO2_concentration = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
	std::vector<double> C7H16_concentration = { 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0 };

	for (auto i = 0; i < CO2_concentration.size(); i++) {
		auto problem = EquilibriumProblem(system);
		problem.setTemperature(394.15);
		problem.setPressure(pressure_range[0]);

		out << std::setw(30) << "T" << std::setw(30) << "P" << std::setw(30) << "converged" << std::setw(30) << "num_phases";

		for (const auto& phase : system.phases())
			out << std::setw(30) << phase.name() << " amount";

		for (const auto& phase : system.phases())
			out << std::setw(30) << phase.name() << " density";

		for (auto i = 0; i < system.species().size(); i++) {
			out << std::setw(15) << system.species(i).name();
		}
		out << std::endl;

		//problem.add("CO2(g)", 1, "mol");
		problem.add("CO2(g)", CO2_concentration[i], "mol");
		problem.add("CH4(g)", C7H16_concentration[i], "mol");

		problem.setPressure(pressure_range[0]);
		auto res = solver.solve(state, problem, phaseidMethod,0);

		auto num_phases_old = count_num_of_phase(state);

		out << "concentration" << std::endl;
		out << "CO2: " << CO2_concentration[i] << std::endl;
		out << "CH4: " << C7H16_concentration[i] << std::endl;
		auto j = 0;
		while (count < 2 && (pressure_range[0] + j * dP) < pressure_range[1]) {
			problem.setPressure(pressure_range[0] + j * dP);
			auto res = solver.solve(state, problem, phaseidMethod,0);
			auto num_phases = count_num_of_phase(state);
			if (num_phases != num_phases_old && res.optimum.succeeded == true) {
				count++;
				output(state, res, out);
				num_phases_old = num_phases;
			}
			j++;
		}
	}





}

void run_alghafri2014_bubble_and_dew_point_CO2_C7H16() {

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
	//editor.addAqueousPhase({ "H2O(l)", "CO2(aq)", "H2S(aq)", "Methane(aq)" });


	// Approach 2:
	std::vector<Phase> phases;
	//phases.push_back(convertPhase(editor.aqueousPhase(), db, pressures, temperatures));

	{
		auto gas_species = std::vector<GaseousSpecies>{
			//db.gaseousSpecies("C2H4(g)"),
			db.gaseousSpecies("CO2(g)"),
			//db.gaseousSpecies("H2O(g)"),
			//db.gaseousSpecies("H2S(g)"),
			//db.gaseousSpecies("O2(g)"),
			db.gaseousSpecies("CH4(g)"),
		};
		auto mixture = GaseousMixture(gas_species);
		auto gas = GaseousPhase(mixture);
		//gas.setChemicalModelRedlichKwong();

		phases.push_back(convertPhase(gas, db));
	}

	{
		auto oil_species = std::vector<HydrocarbonSpecies>{
			HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")),
			//OilSpecies(db.gaseousSpecies("H2S(oil)")),
			//OilSpecies(db.gaseousSpecies("H2O(oil)")),
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
		auto oil = OilPhase(mixture);

		//oil.setChemicalModelRedlichKwong();

		phases.push_back(convertPhase(oil, db));
	}


	// auto system = ChemicalSystem(editor);
	auto system = ChemicalSystem(phases);


	//auto problem = EquilibriumProblem(system);


	phaseBehavior_given_temp(system, { 1e4,100e6 }, "W:\\release\\Projects\\reaktoronote\\notebook\\oilphase\\Z_selection\\phase_diagram\\alghafri_1\\convergence_Z_workanalisys_PR_method_heptane_CO2.txt");

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
    //editor.addAqueousPhase("H O C S");
	editor.addAqueousPhase({  "Methane(aq)" , "CO2(aq)", "H2S(aq)", "H2O(l)" });
	//editor.addAqueousPhase({ "H2O(l)" , "CO2(aq)", "H2S(aq)" });
	//editor.addAqueousPhase({ "H2O(l)"  });
	//editor.addGaseousPhase({"H2O(g)" , "CO2(g)",  "H2S(g)" });// .setChemicalModelRedlichKwong();
	
    //// Approach 1:
    //// editor.addGaseousPhase({"H2O(g)", "CH4(g)" /* , "C2H6(g) " */ });
    //editor.addGaseousPhase({"CH4(g)" /* , "C2H6(g) " */ });
    //ChemicalSystem system(editor);
	
	
    // Approach 2:
    std::vector<Phase> phases;
    phases.push_back(convertPhase(editor.aqueousPhase(), db));

    {
        auto gas_species = std::vector<GaseousSpecies>{
            //db.gaseousSpecies("C2H4(g)"),
			db.gaseousSpecies("CH4(g)"),
			db.gaseousSpecies("CO2(g)"),
			db.gaseousSpecies("H2S(g)"),
			db.gaseousSpecies("H2O(g)"),
			//db.gaseousSpecies("O2(g)"),
			
        };
        auto mixture = GaseousMixture(gas_species);
        auto gas = GaseousPhase(mixture);
		gas.setChemicalModelRedlichKwong();

        phases.push_back(convertPhase(gas, db));
    }

    {
        auto oil_species = std::vector<HydrocarbonSpecies>{
			HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")),
			HydrocarbonSpecies(db.gaseousSpecies("H2O(oil)")),
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
        auto oil = OilPhase(mixture);

		oil.setChemicalModelRedlichKwong();

        phases.push_back(convertPhase(oil, db));
    }


	// auto system = ChemicalSystem(editor);
    auto system = ChemicalSystem(phases);



	//auto p = def_points({ {0, 1.0},{150,120} }, 50);
	//auto p = def_points({ {0, 52.1},{75,120} }, 10);
	auto p = def_points({ {0.0, 1.0} , {210, 240 } }, 100);
	//auto p = std::vector<Point>{ {10, 240} };

	/*
	for (auto i = 0; i < 6333; i++) {
		p.erase(p.begin());
	}
	*/
	std::ofstream out("W:\\release\\Projects\\reaktoronote\\notebook\\\oilphase\\huamg_et_all_validation\\huang_et_all_1985_no_initial_guess_teste.txt");
	
	out << std::setw(30) << "T" << std::setw(30) << "P" << std::setw(30) << "converged" << std::setw(30) << "num_phases";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_amount";

	for (const auto& phase : system.phases())
		out << std::setw(30) << phase.name() << "_density";

	for (auto i = 0; i < system.species().size(); i++) {
		out << std::setw(15) << system.species(i).name();
	}
	out << std::endl;

	auto state = ChemicalState(system);

	for (auto i = 0; i < p.size(); i++) {
		//auto problem = EquilibriumProblem(system);
				
		//const auto& system = problem.system();
		auto system = ChemicalSystem(phases);

		auto problem = EquilibriumProblem(system);

		//problem.setTemperature(p[i].temperature, "degC");
		//problem.setPressure(p[i].pressure, "bar");
		problem.setTemperature(37.8, "degC");
		problem.setPressure(6.26, "MPa");
		problem.add("CO2(g)", 0.0503, "mol");
		problem.add("H2S(g)", 0.3986, "mol");
		problem.add("H2O(g)", 0.5008, "mol");
		//problem.add("O2", 1, "umol");
		//problem.add("H2O", 1, "kg");
		problem.add("CH4(g)", 0.0504, "mol");
		//problem.add("H2S(oil)", 1, "mol");

		auto solver = EquilibriumSolver(problem.system());
		solver.setPartition(problem.partition());

		auto options = EquilibriumOptions();
		options.hessian = GibbsHessian::Exact;
		options.nonlinear.max_iterations = 1000;
		options.optimum.max_iterations = 5000;
		options.optimum.output.active = false;
		options.optimum.ipnewton.step = Conservative;
		options.optimum.tolerance = 1.e-7;
		//options.optimum.tolerance = 1.e-12;

		solver.setOptions(options);

		//if (p[i].pressure == 1.0) {
		state.setSpeciesAmounts(0.0);
		//}
		/*
		state.setSpeciesAmount("H2O(l)", 0.48527);
		state.setSpeciesAmount("CO2(aq)", 0.0017717);
		state.setSpeciesAmount("H2S(aq)", 0.0133178);
		state.setSpeciesAmount("Methane(aq)", 0.000491353);
		state.setSpeciesAmount("CO2(g)", 0.0482283);
		state.setSpeciesAmount("H2O(g)", 0.0145731);
		state.setSpeciesAmount("H2S(g)", 0.386682);
		*/
		//if ((p[i].temperature != p[i - 1].temperature) && i != 0 ) {
		//	state.setSpeciesAmounts(0.0);
		//}
		
		auto phaseidMethod = phaseIdentificationMethod::Gibbs_residual_based;
		auto res = solver.solve(state, problem, phaseidMethod, 0);

		/*
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::VolumeMethod;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod,0);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::workanalisysPengRobinson;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod,0);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::CriticalPointMethods;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod,0);
		}
		if (res.optimum.succeeded == false) {
			phaseidMethod = phaseIdentificationMethod::workanalisysPengRobinson;
			state.setSpeciesAmounts(0.0);
			res = solver.solve(state, problem, phaseidMethod,0);
		}
		*/

		//std::cout << state << std::endl;

		//auto final_res = solver.solve(state, problem, phaseidMethod, -190);
		/*
		if (res.optimum.succeeded == false) {
			state.setSpeciesAmounts(0.0);
			auto res = solver.solve(state, problem, phaseidMethod, 0);
		}
		*/
		output(state, res, out);

		//if (res.optimum.succeeded == false) {
		//	state.setSpeciesAmounts(0.0);
		//}
		std::cout << state << std::endl;
		//std::cout << final_res.optimum.succeeded << std::endl;
	}

	/*
	auto state = ChemicalState(system);

	auto solver = EquilibriumSolver(system);		

	auto options = EquilibriumOptions();
	options.hessian = GibbsHessian::Exact;
	options.nonlinear.max_iterations = 100;
	options.optimum.max_iterations = 200;

	solver.setOptions(options);
	
	//state.setSpeciesAmounts(0.0);

	auto outputoption = OutputterOptions();

	outputoption = true;

	auto outputter = Outputter();
	
	auto phaseidMethod = phaseIdentificationMethod::Gibbs_residual_based;

	outputter.setOptions(outputoption);
	auto res = solver.solve(state, problem, phaseidMethod);
	std::cout << state << std::endl;

	std::cout << res.optimum.succeeded << std::endl;
	
	phaseidMethod = phaseIdentificationMethod::Gibbs_residual_based;
	res = solver.solve(state, problem, phaseidMethod);
	
	

	std::cout << state << std::endl;

	std::cout << res.optimum.succeeded << std::endl;
	
	
    */
	//auto p = def_points({ {10, 1},{160,400} }, 50);
	
	//distance_order(p);

	//build_flash(problem, p, "W:\\release\\Projects\\reaktoronote\\notebook\\oilphase\\Z_selection\\convergence_study\\oil_gas_lump_CO2_H2S\\convergence_Z_workanalisys_PR_method_oil_gas_lump_CO2_H2S.txt");

	//phaseBehavior_given_temp(problem, { 0, 15e6 }, "W:\\release\\Projects\\reaktoronote\\notebook\\oilphase\\Z_selection\\convergence_study\\mixture4\\convergence_Z_workanalisys_PR_method_oil_gas_lump_CO2_H2S.txt");

	//run_mixture_Huang_1984();

	//run_alghafri2014_bubble_and_dew_point_CO2_C7H16();

}