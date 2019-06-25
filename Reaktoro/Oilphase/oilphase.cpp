#include <Reaktoro/Reaktoro.hpp>

using namespace Reaktoro;


// -------------------------------------------------------------------------------------------------

// NOTE: Copied (and modified) from GaseousPhase, GaseousMixture, GaseousSpecies, etc.:


/// A type to describe the attributes of a gaseous species
class OilSpecies : public Species
{
private:
    struct Impl
    {
        // The critical temperature of the gaseous species (in units of K)
        double critical_temperature = 0.0;

        // The critical pressure of the gaseous species (in units of Pa)
        double critical_pressure = 0.0;

        // The acentric factor of the gaseous species
        double acentric_factor = 0.0;

        /// The thermodynamic data of the gaseous species.
        HydrocarbonSpeciesThermoData thermo;
    };

public:
    /// Construct a default OilSpecies instance
    OilSpecies();

    /// Construct an OilSpecies instance from a Species instance
    OilSpecies(const Species& species);

    // [Prototyping: added now]
    /// Copies data from a GaseousSpecies
    explicit OilSpecies(const GaseousSpecies& species)
        : Species(species)
        , pimpl(new Impl())
    {
        this->setCriticalTemperature(species.criticalTemperature());
        this->setCriticalPressure(species.criticalPressure());
        this->setAcentricFactor(species.acentricFactor());

        auto const& gas_thermo_data = species.thermoData();
        
        Optional<HydrocarbonSpeciesThermoParamsHKF> oil_hkf;
        if (!gas_thermo_data.hkf.empty()) {
            GaseousSpeciesThermoParamsHKF const& hkf = gas_thermo_data.hkf.get();
            oil_hkf = HydrocarbonSpeciesThermoParamsHKF { hkf.Gf, hkf.Hf, hkf.Sr, hkf.a, hkf.b, hkf.c, hkf.Tmax };
        }

        this->setThermoData(HydrocarbonSpeciesThermoData{
            gas_thermo_data.properties,
            gas_thermo_data.reaction,
            oil_hkf,
            gas_thermo_data.phreeqc
        });
    }

    // [End]

    /// Set the critical temperature of the gaseous species (in units of K)
    auto setCriticalTemperature(double val) -> void;

    /// Set the critical pressure of the gaseous species (in units of Pa)
    auto setCriticalPressure(double val) -> void;

    /// Set the acentric factor of the gaseous species
    auto setAcentricFactor(double val) -> void;

    /// Set the thermodynamic data of the gaseous species.
    auto setThermoData(const HydrocarbonSpeciesThermoData& thermo) -> void;

    /// Return the critical temperature of the gaseous species (in units of K)
    auto criticalTemperature() const -> double;

    /// Return the critical pressure of the gaseous species (in units of Pa)
    auto criticalPressure() const -> double;

    /// Return the acentric factor of the gaseous species
    auto acentricFactor() const -> double;

    /// Return the thermodynamic data of the gaseous species.
    auto thermoData() const -> const HydrocarbonSpeciesThermoData&;

private:
    std::shared_ptr<Impl> pimpl;
};


OilSpecies::OilSpecies()
: pimpl(new Impl())
{}

OilSpecies::OilSpecies(const Species& species)
: Species(species), pimpl(new Impl())
{}

auto OilSpecies::setCriticalTemperature(double val) -> void
{
    Assert(val > 0.0, "Cannot set the critical temperature of the gas `" + name() + "`.",
        "The given critical temperature `" + std::to_string(val) + "` is not positive.");
    pimpl->critical_temperature = val;
}

auto OilSpecies::setCriticalPressure(double val) -> void
{
    Assert(val > 0.0, "Cannot set the critical pressure of the gas `" + name() + "`.",
        "The given critical pressure `" + std::to_string(val) + "` is not positive.");
    pimpl->critical_pressure = val;
}

auto OilSpecies::setAcentricFactor(double val) -> void
{
    pimpl->acentric_factor = val;
}

auto OilSpecies::setThermoData(const HydrocarbonSpeciesThermoData& thermo) -> void
{
    pimpl->thermo = thermo;
}

auto OilSpecies::criticalTemperature() const -> double
{
    return pimpl->critical_temperature;
}

auto OilSpecies::criticalPressure() const -> double
{
    return pimpl->critical_pressure;
}

auto OilSpecies::acentricFactor() const -> double
{
    return pimpl->acentric_factor;
}

auto OilSpecies::thermoData() const -> const HydrocarbonSpeciesThermoData&
{
    return pimpl->thermo;
}


/// Provides a computational representation of a gaseous mixture.
/// The OilMixture class is defined as a collection of OilSpecies objects,
/// representing, therefore, a mixture of gaseous species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of gaseous
/// species.
/// @see OilSpecies
/// @ingroup Mixtures
class OilMixture : public GeneralMixture<HydrocarbonSpecies>
{
public:
    /// Construct a default OilMixture instance.
    OilMixture();

    /// Construct a OilMixture instance with given species.
    /// @param species The species that compose the gaseous mixture
    explicit OilMixture(const std::vector<HydrocarbonSpecies>& species);

    /// Destroy the OilMixture instance.
    virtual ~OilMixture();

    /// Calculate the state of the gaseous mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, VectorConstRef n) const -> MixtureState;
};

OilMixture::OilMixture()
: GeneralMixture<HydrocarbonSpecies>()
{}

OilMixture::OilMixture(const std::vector<HydrocarbonSpecies>& species)
: GeneralMixture<HydrocarbonSpecies>(species)
{}

OilMixture::~OilMixture()
{}

auto OilMixture::state(Temperature T, Pressure P, VectorConstRef n) const -> MixtureState
{
    MixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}


/// Class that defines a gaseous phase
class OilPhase : public Phase
{
public:
    /// Construct a default OilPhase instance.
    OilPhase();

    /// Construct an OilPhase instance with given gaseous mixture.
    /// The Peng-Robinson equation of state is chosen by default to calculate the
    /// thermodynamic and chemical properties of this OilPhase object.
    explicit OilPhase(const OilMixture& mixture);

    /// Set the chemical model of the phase with the ideal gas equation of state.
    auto setChemicalModelIdeal() -> OilPhase&;

    /// Set the chemical model of the phase with the van der Waals equation of state.
    /// Reference: *van der Waals, J.D. (1910). The equation of state for gases and liquids. Nobel Lectures in Physics. pp. 254-265*.
    auto setChemicalModelVanDerWaals() -> OilPhase&;

    /// Set the chemical model of the phase with the Redlich-Kwong equation of state.
    /// Reference: *Redlich, O., Kwong, J.N.S. (1949). On The Thermodynamics of Solutions. Chem. Rev. 44(1) 233–244*.
    auto setChemicalModelRedlichKwong() -> OilPhase&;

    /// Set the chemical model of the phase with the Soave-Redlich-Kwong equation of state.
    /// Reference: *Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci., 27, 1197-1203*.
	auto setChemicalModelSoaveRedlichKwong()->OilPhase&;

    /// Set the chemical model of the phase with the Peng-Robinson equation of state.
    /// Reference: *Peng, D.Y., Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial and Engineering Chemistry: Fundamentals 15: 59–64*.
    auto setChemicalModelPengRobinson() -> OilPhase&;

    /// Set the chemical model of the phase with the Spycher et al. (2003) equation of state.
    /// This model only supports the gaseous species `H2O(g)` and `CO2(g)`. Any other species
    /// will result in a runtime error.
    /// Reference: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
    /// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
    /// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
    auto setChemicalModelSpycherPruessEnnis() -> OilPhase&;

    /// Set the chemical model of the phase with the Spycher and Reed (1988) equation of state.
    /// This model only supports the gaseous species `H2O(g)`, `CO2(g)`, and CH4(g). Any other
    /// species will result in a runtime error.
    /// Reference: *Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2,
    /// CH4, H2O and of H2O--CO2--CH4 mixtures: A virial equation treatment for
    /// moderate pressures and temperatures applicable to calculations of
    /// hydrothermal boiling. Geochimica et Cosmochimica Acta, 52(3), 739–749*.
    auto setChemicalModelSpycherReed() -> OilPhase&;

    /// Return the OilMixture instance
    auto mixture() const -> const OilMixture&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};


struct OilPhase::Impl
{
    /// The gaseous mixture instance
    OilMixture mixture;

    /// Construct a default Impl instance
    Impl()
    {
    }

    /// Construct a custom Impl instance
    Impl(const OilMixture& mixture)
        : mixture(mixture)
    {
    }
};

OilPhase::OilPhase()
    : Phase(), pimpl(new Impl())
{
}

OilPhase::OilPhase(const OilMixture& mixture)
    : pimpl(new Impl(mixture))
{
    // Convert the HydrocarbonSpecies instances to Species instances
    std::vector<Species> species;
    for (const HydrocarbonSpecies& x : mixture.species())
        species.push_back(x);

    // Set the Phase attributes
    setName("Oil");
    setType(PhaseType::Gas);
    setSpecies(species);
    setChemicalModelPengRobinson();
}

//auto OilPhase::setChemicalModelIdeal() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelIdeal(mixture());
//    setChemicalModel(model);
//    return *this;
//}
//
//auto OilPhase::setChemicalModelVanDerWaals() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelVanDerWaals(mixture());
//    setChemicalModel(model);
//    return *this;
//}
//
//auto OilPhase::setChemicalModelRedlichKwong() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelRedlichKwong(mixture());
//    setChemicalModel(model);
//    return *this;
//}
//
//auto OilPhase::setChemicalModelSoaveRedlichKwong() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelSoaveRedlichKwong(mixture());
//    setChemicalModel(model);
//    return *this;
//}
//

auto oilChemicalModelCubicEOS(const OilMixture& mixture, CubicEOS::Model modeltype)->PhaseChemicalModel
{
    // Copy & Paste
    // The number of gases in the mixture
    const unsigned nspecies = mixture.numSpecies();

    // Get the the critical temperatures, pressures and acentric factors of the gases
    std::vector<double> Tc, Pc, omega;
    for(HydrocarbonSpecies species : mixture.species())
    {
        Tc.push_back(species.criticalTemperature());
        Pc.push_back(species.criticalPressure());
        omega.push_back(species.acentricFactor());
    }

    // Initialize the CubicEOS instance
    CubicEOS eos(nspecies);

    eos.setPhaseAsLiquid();  // <-
    // eos.setPhaseAsVapor();  // <-
	//eos.setMinimazeGibbsEnergytrue();

    eos.setCriticalTemperatures(Tc);
    eos.setCriticalPressures(Pc);
    eos.setAcentricFactors(omega);
    eos.setModel(modeltype);

    // The state of the mixture
    MixtureState state;

    // Define the chemical model function of the gaseous phase
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // Evaluate the state of the gaseous mixture
        state = mixture.state(T, P, n);

        // The mole fractions of the species
        const auto& x = state.x;

        // Evaluate the CubicEOS object function
        const CubicEOS::Result eosres = eos(T, P, x);

        // The ln of mole fractions
        const ChemicalVector ln_x = log(x);

        // The ln of pressure in bar units
        const ThermoScalar ln_Pbar = log(1e-5 * P);

        // Create an alias to the ln fugacity coefficients
        const auto& ln_phi = eosres.ln_fugacity_coefficients;

        // Fill the chemical properties of the gaseous phase
        res.ln_activity_coefficients = ln_phi;
        res.ln_activities = ln_phi + ln_x + ln_Pbar;
        res.molar_volume = eosres.molar_volume;
        res.residual_molar_gibbs_energy = eosres.residual_molar_gibbs_energy;
        res.residual_molar_enthalpy = eosres.residual_molar_enthalpy;
        res.residual_molar_heat_capacity_cp = eosres.residual_molar_heat_capacity_cp;
        res.residual_molar_heat_capacity_cv = eosres.residual_molar_heat_capacity_cv;

		/*
		//make activities of phase that did not supoused to be ther equal to zero		
		bool should_be_zero = true;
		for (auto i = 0; i < res.ln_activity_coefficients.val.size(); i++) {
			if (res.ln_activity_coefficients.val[i] != 0.0) {
				should_be_zero = false;
			}
		}
		if (should_be_zero) {
			res.ln_activities.fill(0.0);
		}
		*/
		/*
		for (int i = 0; i < res.ln_activities.val.size(); i++) {
			for (int j = 0; j < res.ln_activities.val.size(); j++) {
				std::cout << res.ln_activities.ddn(i, j) << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;
		*/
    };

    return model;
}



auto oilChemicalModelPengRobinson(const OilMixture& mixture) -> PhaseChemicalModel
{
    return oilChemicalModelCubicEOS(mixture, CubicEOS::PengRobinson);
}



auto OilPhase::setChemicalModelPengRobinson() -> OilPhase&
{
    PhaseChemicalModel model = oilChemicalModelPengRobinson(mixture());
    setChemicalModel(model);
    return *this;
}

auto oilChemicalModelRedlichKwong(const OilMixture& mixture) -> PhaseChemicalModel
{
	return oilChemicalModelCubicEOS(mixture, CubicEOS::RedlichKwong);
}



auto OilPhase::setChemicalModelRedlichKwong() -> OilPhase&
{
	PhaseChemicalModel model = oilChemicalModelRedlichKwong(mixture());
	setChemicalModel(model);
	return *this;
}

auto oilChemicalModelSoaveRedlichKwong(const OilMixture& mixture) -> PhaseChemicalModel
{
	return oilChemicalModelCubicEOS(mixture, CubicEOS::SoaveRedlichKwong);
}



auto OilPhase::setChemicalModelSoaveRedlichKwong() -> OilPhase&
{
	PhaseChemicalModel model = oilChemicalModelSoaveRedlichKwong(mixture());
	setChemicalModel(model);
	return *this;
}


//
//auto OilPhase::setChemicalModelSpycherPruessEnnis() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelSpycherPruessEnnis(mixture());
//    setChemicalModel(model);
//    return *this;
//}
//
//auto OilPhase::setChemicalModelSpycherReed() -> OilPhase&
//{
//    PhaseChemicalModel model = gaseousChemicalModelSpycherReed(mixture());
//    setChemicalModel(model);
//    return *this;
//}

auto OilPhase::mixture() const -> const OilMixture&
{
    return pimpl->mixture;
}



// -------------------------------------------------------------------------------------------------




// NOTE: Copied from ChemicalEditor [PROTOTYPING]

auto lnActivityConstants(const AqueousPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the aqueous species
    ThermoVector ln_c(phase.numSpecies());

    // The index of solvent water species
    const Index iH2O = phase.indexSpeciesAnyWithError(alternativeWaterNames());

    // Set the ln activity constants of aqueous species to ln(55.508472)
    ln_c = std::log(1.0 / waterMolarMass);

    // Set the ln activity constant of water to zero
    ln_c[iH2O] = 0.0;

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        return ln_c;
    };

    return f;
}

auto lnActivityConstants(const GaseousPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the gaseous species
    ThermoVector ln_c(phase.numSpecies());

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        ln_c = log(P * 1e-5); // ln(Pbar)
        return ln_c;
    };

    return f;
}

auto lnActivityConstants(const OilPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the gaseous species
    ThermoVector ln_c(phase.numSpecies());

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        ln_c = log(P * 1e-5); // ln(Pbar)
        return ln_c;
    };

    return f;
}



// NOTE: Copy & Pasted from ChemicalEditor [PROTOTYPING]
template<typename PhaseType>
auto convertPhase(const PhaseType& phase, const Database& database) -> Phase
{
	// NOTE: Copy & Pasted from ChemicalEditor {
	// The default temperatures for the interpolation of the thermodynamic properties (in units of celsius)
	std::vector<double> temperatures{ 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300 };

	// The default pressures for the interpolation of the thermodynamic properties (in units of bar)
	std::vector<double> pressures{ 1, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 };

	// Convert the temperatures and pressures to units of kelvin and pascal respectively
	for (auto& x : temperatures) { x = x + 273.15; }
	for (auto& x : pressures) { x = x * 1.0e+5; }

	// The number of species in the phase
	const unsigned nspecies = phase.numSpecies();

	// Define the lambda functions for the calculation of the essential thermodynamic properties
	Thermo thermo(database);

	std::vector<ThermoScalarFunction> standard_gibbs_energy_fns(nspecies);
	std::vector<ThermoScalarFunction> standard_enthalpy_fns(nspecies);
	std::vector<ThermoScalarFunction> standard_volume_fns(nspecies);
	std::vector<ThermoScalarFunction> standard_heat_capacity_cp_fns(nspecies);
	std::vector<ThermoScalarFunction> standard_heat_capacity_cv_fns(nspecies);

	// Create the ThermoScalarFunction instances for each thermodynamic properties of each species
	for (unsigned i = 0; i < nspecies; ++i)
	{
		const std::string name = phase.species(i).name();

		standard_gibbs_energy_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarGibbsEnergy(T, P, name); };
		standard_enthalpy_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarEnthalpy(T, P, name); };
		standard_volume_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarVolume(T, P, name); };
		standard_heat_capacity_cp_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstP(T, P, name); };
		standard_heat_capacity_cv_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstV(T, P, name); };
	}

	// Create the interpolation functions for thermodynamic properties of the species
	ThermoVectorFunction standard_gibbs_energies_interp = interpolate(temperatures, pressures, standard_gibbs_energy_fns);
	ThermoVectorFunction standard_enthalpies_interp = interpolate(temperatures, pressures, standard_enthalpy_fns);
	ThermoVectorFunction standard_volumes_interp = interpolate(temperatures, pressures, standard_volume_fns);
	ThermoVectorFunction standard_heat_capacities_cp_interp = interpolate(temperatures, pressures, standard_heat_capacity_cp_fns);
	ThermoVectorFunction standard_heat_capacities_cv_interp = interpolate(temperatures, pressures, standard_heat_capacity_cv_fns);
	ThermoVectorFunction ln_activity_constants_func = lnActivityConstants(phase);

	// Define the thermodynamic model function of the species
	PhaseThermoModel thermo_model = [=](PhaseThermoModelResult& res, Temperature T, Pressure P)
	{
		// Calculate the standard thermodynamic properties of each species
		res.standard_partial_molar_gibbs_energies = standard_gibbs_energies_interp(T, P);
		res.standard_partial_molar_enthalpies = standard_enthalpies_interp(T, P);
		res.standard_partial_molar_volumes = standard_volumes_interp(T, P);
		res.standard_partial_molar_heat_capacities_cp = standard_heat_capacities_cp_interp(T, P);
		res.standard_partial_molar_heat_capacities_cv = standard_heat_capacities_cv_interp(T, P);
		res.ln_activity_constants = ln_activity_constants_func(T, P);

		return res;
	};

	// Create the Phase instance
	Phase converted = phase;
	converted.setThermoModel(thermo_model);

	return converted;
}
