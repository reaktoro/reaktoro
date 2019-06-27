#include <Reaktoro/Reaktoro.hpp>

#include <math.h>

using namespace Reaktoro;



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
	phases.push_back(editor.convertAqueousPhase(editor.aqueousPhase()));

	{
		auto gas_species = std::vector<GaseousSpecies>{
			db.gaseousSpecies("CO2(g)"),
		};
		auto mixture = GaseousMixture(gas_species);
		auto gas = GaseousPhase(mixture);

		phases.push_back(editor.convertGaseousPhase(gas));
	}

	{
		auto oil_species = std::vector<HydrocarbonSpecies>{
			HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")),
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

		//oil.setChemicalModelSoaveRedlichKwong();

		//oil.setChemicalModelRedlichKwong();

		phases.push_back(editor.convertHydrocarbonPhase(oil));
	}


	
	auto system = ChemicalSystem(phases);

	auto problem = EquilibriumProblem(system);

	problem.setTemperature(-10 , "degC");
	problem.setPressure(20, "bar");
	problem.add("CH4(g)", 1, "mol");

	auto gas_system = ChemicalSystem({system.phase("Gaseous")});
	auto liq_system = ChemicalSystem({system.phase("Oil")});

	auto z = Eigen::VectorXd(1);
	auto K_liq = Eigen::VectorXd(1);
	auto K_vap = Eigen::VectorXd(1);
	auto Y_liq = Eigen::VectorXd(1);
	auto y_liq = Eigen::VectorXd(1);
	auto Y_vap = Eigen::VectorXd(1);
	auto y_vap = Eigen::VectorXd(1);

	auto R_liq = Eigen::VectorXd(1);
	auto R_vap = Eigen::VectorXd(1);

	z(0) = 1.0;
	
	auto z_cres_gas = ChemicalModelResult(1, 1);
	auto z_cres_liq = ChemicalModelResult(1, 1);
	auto z_cres = ChemicalModelResult(1, 1);
	gas_system.chemicalModel()(z_cres_gas, problem.temperature(), problem.pressure(), z);
	liq_system.chemicalModel()(z_cres_liq, problem.temperature(), problem.pressure(), z);

	//if (z_cres_gas.phaseResidualMolarGibbsEnergies()[0].val > z_cres_liq.phaseResidualMolarGibbsEnergies()[0].val) {
	z_cres = z_cres_liq;
	//}

	
	for (auto i = 0; i < gas_system.numSpecies(); i++) {
		const auto& Pc = db.gaseousSpecies(gas_system.species(i).name()).criticalPressure();// .critical_temperature();
		const auto& Tc = db.gaseousSpecies(gas_system.species(i).name()).criticalTemperature();
		const auto& w = db.gaseousSpecies(gas_system.species(i).name()).acentricFactor();
		auto Ki = (Pc / problem.pressure())*std::exp(5.37 * (1.0 + w)*(1.0 - (Tc/problem.temperature())));
		K_liq(i) = Ki;
		K_vap(i) = Ki;
	}
	
	bool converged = false;
	bool TSliq = false;
	bool TSvap = false;

	double Sl = 0.0;
	double Sv = 0.0;
	for (auto i = 0; i < 10000; i++) {
		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			Y_liq(j) = z(j)/K_liq(j);
			Y_vap(j) = z(j)*K_vap(j);
		}
		
		Sl = Y_liq.sum();
		Sv = Y_vap.sum();

		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			y_liq(j) = Y_liq(j) / Sl;
			y_vap(j) = Y_vap(j) / Sv;
		}

		auto y_liq_cres = ChemicalModelResult(1, 1);
		auto y_vap_cres = ChemicalModelResult(1, 1);

		liq_system.chemicalModel()(y_liq_cres, problem.temperature(), problem.pressure(), y_vap);
		gas_system.chemicalModel()(y_vap_cres, problem.temperature(), problem.pressure(), y_liq);

		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			R_liq(j) =  (Sl)*(std::exp(y_liq_cres.lnActivityCoefficients()[j].val)/std::exp(z_cres.lnActivityCoefficients()[j].val));
			R_vap(j) = (1.0/Sv)*(std::exp(z_cres.lnActivityCoefficients()[j].val)/std::exp(y_vap_cres.lnActivityCoefficients()[j].val));
		}

		auto trial_solution_liq = 0.0;
		auto trial_solution_vap = 0.0;
		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			trial_solution_liq += std::log(K_liq(j))*std::log(K_liq(j));
			trial_solution_vap += std::log(K_vap(j))*std::log(K_vap(j));
		}
		if (trial_solution_liq < 1.e-4) {
			TSliq = true;
		}
		if (trial_solution_vap < 1.e-4) {
			TSvap = true;
		}

		auto erro_liq = 0.0;
		auto erro_vap = 0.0;
		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			erro_liq += (R_liq(j) - 1)*(R_liq(j) - 1);
			erro_vap += (R_vap(j) - 1)*(R_vap(j) - 1);
		}
		if (erro_liq < 1e-12) {
			std::cout << "liq possivel to form" << std::endl;
			converged = true;
		}
		if (erro_vap < 1e-12) {
			std::cout << "vap to form" << std::endl;
			converged = true;
		}
		if (converged) {
			break;
		}

		for (auto j = 0; j < gas_system.numSpecies(); j++) {
			K_liq(j) = K_liq(j) * R_liq(j);
			K_vap(j) = K_vap(j) * R_liq(j);
		}
		if (i == 99999) {
			std::cout << "max points" << std::endl;
		}
	}

	if (converged){
		if (TSliq && TSvap) {
			std::cout << "1 fase" << std::endl;
		}else if (Sl <= 1.0 && TSvap) {
			std::cout << "2 fases" << std::endl;
		}else if (TSliq && Sl <= 1.0) {
			std::cout << "2 fases" << std::endl;
		}else if (Sv <= 1.0 && Sl <= 1.0) {
			std::cout << "3 fases" << std::endl;
		}
	}else {
		if (Sv > 1 && TSliq) {
			std::cout << "2 fases" << std::endl;
		} else if (TSvap && Sl > 1) {
			std::cout << "2 fases" << std::endl;
		} else if (Sv > 1 && Sl > 1) {
			std::cout << "2 fases" << std::endl;
		} else if (Sv > 1 && Sl <= 1) {
			std::cout << "3 fases" << std::endl;
		} else if (Sv <= 1 && Sl > 1) {
			std::cout << "3 fases" << std::endl;
		}
	}

	std::cout << K_liq(0) << std::endl;
	std::cout << K_vap(0) << std::endl;
	std::cout << "ok" << std::endl;
	
	

}