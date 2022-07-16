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

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

// Reactive transport test includes
#include <demos/cpp/DemosUtils.h>

using namespace Reaktoro;

/// Forward declaration
auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportKineticsResults& results) -> void;

int main()
{
    tic(TOTAL_TIME)

    // Step 1: Initialise auxiliary time-related constants
    // int second = 1;
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;
    int week = 7 * day;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 2000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 60.0;                     // the temperature (in units of degC)
    params.P = 100;                      // the pressure (in units of bar)

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // ----------------------------------------------------------- //
    // Smart equilibrium parameters
    // ----------------------------------------------------------- //

    // Run smart algorithm with clustering
    params.method = SmartEquilibriumStrategy::Clustering;
    params.smart_equilibrium_reltol = 1e-3;

    // ----------------------------------------------------------- //
    // Smart kinetic parameters
    // ----------------------------------------------------------- //

    // Select clustering approach
    params.smart_kin_method = SmartKineticStrategy::Clustering;
    params.smart_kinetic_tol = 1e-3; // 5e-3; // 1e-4;
    params.smart_kinetic_reltol = 1e-1;
    params.smart_kinetic_abstol = 1e-4;

//    // Select priority queue approach
//    params.smart_kin_method = SmartKineticStrategy::PriorityQueue;
//    params.smart_kinetic_tol = 1e-3;
//    params.smart_kinetic_abstol = 1e-4;
//    params.smart_kinetic_reltol = 1e-1;
//
//    // Select nearest neighbour approach
//    params.smart_kin_method = SmartKineticStrategy::NearestNeighbour;
//    params.smart_kinetic_tol = 1e-3;
//    params.smart_kinetic_abstol = 1e-4;
//    params.smart_kinetic_reltol = 1e-4;


    // ----------------------------------------------------------- //
    // Optimization algorithm parameters
    // ----------------------------------------------------------- //

    params.hessian = GibbsHessian::Exact;
    //params.hessian = GibbsHessian::Approximation;

    // ----------------------------------------------------------- //
    // Activity model for the aqueous species
    // ----------------------------------------------------------- //

    // Define the activity model for the aqueous species
    params.activity_model = ActivityModel::HKF;
    //params.activity_model = ActivityModel::HKFSelectedSpecies;
    //params.activity_model = ActivityModel::Pitzer;
    //params.activity_model = ActivityModel::PitzerSelectedSpecies;
    //params.activity_model = ActivityModel::DebyeHuckel;
    //params.activity_model = ActivityModel::DebyeHuckelSelectedSpecies;

    // Output
    params.outputConsole();

    // Results
    ReactiveTransportKineticsResults results;

    /// **************************************************************************************************************///
    /// CONVENTIONAL kinetics & SMART equilibrium
    /// **************************************************************************************************************///
    /// params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);

    /// **************************************************************************************************************///
    //// SMART kinetics & CONVENTIONAL equilibrium
    /// **************************************************************************************************************///
//    params.use_smart_kinetics_solver = true; params.use_smart_equilibrium_solver = false;
//    params.outputConsoleKineticMethod();
//    runReactiveTransport(params, results);

    /// **************************************************************************************************************///
    /// SMART kinetics & SMART equilibrium
    /// **************************************************************************************************************///
    params.use_smart_kinetics_solver = true;  params.use_smart_equilibrium_solver = true;
    params.outputConsoleKineticMethod();
    runReactiveTransport(params, results);

//    /// **************************************************************************************************************///
//    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
//    /// **************************************************************************************************************///
//    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = false;
//    params.outputConsoleKineticMethod();
//    runReactiveTransport(params, results);

    // **************************************************************************************************************///
    // SPEED-UP analysis
    // **************************************************************************************************************///
    if(results.smart_kin_conv_eq_total != 0){
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart kinetics solver                     : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_conv_eq_total_ideal_search_store_properties << std::endl;

        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT smart.kin.& conv.eq. : " << results.time_reactive_transport_smart_kin_conv_eq << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_smart_kin_conv_eq << std::endl;
    }
    if(results.conv_kin_smart_eq_total != 0) {
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart equilibrium solver                  : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.conv_kin_smart_eq_total_ideal_search_store_properties
                  << std::endl;

        std::cout << "speed up in equilibration    : "
                  << results.conv_kin_conv_eq_total_equilibration /
                     results.conv_kin_smart_eq_total_smart_equilibration << std::endl;
        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT conv.kin.& smart.eq. : " << results.time_reactive_transport_conv_kin_smart_eq
                  << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_conv_kin_smart_eq << std::endl;
    }
    if(results.smart_kin_smart_eq_total != 0){
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "speed up of using smart kinetic solver (smart equilibrium)  : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total << std::endl;
        std::cout << "speed up ... (with ideal search)                            : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search << std::endl;
        std::cout << "speed up ... (with ideal search & store)                    : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search_store << std::endl;
        std::cout << "speed up ... (with ideal search & store & properties eval.) : "
                  << results.conv_kin_conv_eq_total / results.smart_kin_smart_eq_total_ideal_search_store_properties << std::endl;

        std::cout << "time RT conv.kin.& conv.eq.  : " << results.time_reactive_transport_conv_kin_conv_eq << std::endl;
        std::cout << "time RT smart.kin.& smart.eq. : " << results.time_reactive_transport_smart_kin_smart_eq << std::endl;
        std::cout << "speedup                      : " << results.time_reactive_transport_conv_kin_conv_eq
                                                          / results.time_reactive_transport_smart_kin_smart_eq << std::endl;
    }
    std::cout << "total time : " << toc(TOTAL_TIME) << std::endl;

    return 0;
}

auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportKineticsResults& results) -> void
{

    // Step **: Create the results folder
    auto folder = params.makeResultsFolderKinetics("dolomitization-calcite-dolomite");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = params.hessian;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.learning.hessian = params.hessian;
    smart_equilibrium_options.method = params.method;

    // Step **: Define chemical kinetic solver options
    KineticOptions kinetic_options;
    kinetic_options.equilibrium = equilibrium_options;
    kinetic_options.smart_equilibrium = smart_equilibrium_options;
    kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;

    // Step **: Define smart chemical kinetic solver options
    SmartKineticOptions smart_kinetic_options;
    smart_kinetic_options.reltol = params.smart_kinetic_reltol;
    smart_kinetic_options.abstol = params.smart_kinetic_abstol;
    smart_kinetic_options.tol = params.smart_kinetic_tol;
    smart_kinetic_options.learning = kinetic_options;
    smart_kinetic_options.learning.equilibrium = equilibrium_options;
    smart_kinetic_options.smart_equilibrium = smart_equilibrium_options;
    smart_kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    smart_kinetic_options.method = params.smart_kin_method;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;

    // Define the list of selected elements
    StringList selected_elements = "H O Na Cl Ca Mg C";

    // Define the list of selected species
    StringList selected_species = "H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3-- CaCl+ Ca(HCO3)+ MgCl+ Mg(HCO3)+";

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == ActivityModel::HKF){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == ActivityModel::HKFSelectedSpecies){
        // HKF selected species
        editor.addAqueousPhase(selected_species);
    }
    else if(params.activity_model == ActivityModel::Pitzer){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ActivityModel::PitzerSelectedSpecies){
        // Pitzer selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ActivityModel::DebyeHuckel){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ActivityModel::DebyeHuckelSelectedSpecies){
        // Debye-Huckel selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    // Step **: Add mineral phase
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Quartz");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");
    reaction.setSpecificSurfaceArea(50, "cm2/g");

    // Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Ionic strength function
    const auto I = ChemicalProperty::ionicStrength(system);

    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The universal gas constant (in units of kJ/(mol*K))
    const double R = 8.3144621e-3;

    // Time step
    double dt = params.dt;

    //Calcite
    //CaCO3 + 1.000H+ = 1.000HCO3- + 1.000Ca+2
    //     log_k     1.847
    //     delta_h  -25.325    #kJ/mol        #Internal calculation
    //     -analytic -8.5010157E+2  -1.3947146E-1  4.6881027E+4  3.0964897E+2  -2.6591521E+6
    //     #References = LogK/DGf: 06bla/pia; DHf/DHr: Internal calculation; S°: 82plu/bus; Cp: 95rob/hem; V°: 78hel/del,82plu/bus;
    std::string eq_str_calcite = "Calcite + H+  =  Ca++ + HCO3-";
    MineralReaction min_reaction_calcite = editor.addMineralReaction("Calcite")
            .setEquation(eq_str_calcite)
            .setSpecificSurfaceArea(0.7, "m2/g"); // S = 0.7 in m2/g
    Reaction reaction_calcite = createReaction(min_reaction_calcite, system);
    reaction_calcite.setName("Calcite reaction");
    ReactionRateFunction rate_func_calcite_olimse_dat = [&min_reaction_calcite, &reaction_calcite, &system, &num_species, &R, &dt](const ChemicalProperties& properties) -> ChemicalScalar {

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0), res_growth(num_species, 0.0), res_nuc(num_species, 0.0);

        // Auxiliary variables with chemical scalars of value one
        ChemicalScalar f(num_species, 1.0);

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_calcite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_calcite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0 = reaction_calcite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_calcite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        // Calculate  MOL("Ca+2")
        // const Index i_ca = system.indexSpeciesWithError("Ca++");
        // const auto mol_ca = n[i_ca].val;
        // Calculate TOT("HCO3-")
        // const Index i_hco3 = system.indexSpeciesWithError("HCO3-");
        // const auto mol_hco3 = n[i_hco3].val;

        // Calculate activities for H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_hco3 = system.indexSpeciesWithError("HCO3-");
        double activity_hco3 = std::exp(lna(i_hco3));

        // The molar mass of the mineral (in units of kg/mol) # Mm = 100.087
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 350
        // If (SRmin = 1) Then GoTo 350
        if((nm <= 0 && Omega < 1) || Omega == 1) // the is no way to dissolve further
            return res;

        // Average BET
        const auto ssa = min_reaction_calcite.specificSurfaceArea();

        // Molar volume in m3/mol # Mv = 3.693E-05
        const auto molar_volume = 3.693E-05;

        if(Omega < 1) // Omega < 1, dissolution kinetics # kinetic data extracted from 04pal/kha
        {
            // Neutral mechanism # knu = 1.6E-6 * exp((-24000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 1.6E-6 * exp(- 24.0 / R * (1.0 / T - 1.0 / 298.15));

            // Acidic mechanism # k1 = 5E-1 * exp((-14000 / 8.314) * ((1 / TK) - (1 / 298.15))) * ACT("H3O+")
            const auto kappa_acid = 5E-1 * exp(- 14.0 / R * (1.0 / T - 1.0 / 298.15)) * activity_h;

            // total rate constant
            const auto kappa_diss = kappa_neu + kappa_acid;

            // Calculate the resulting mechanism function
            // rate = (knu + k1) * S * m * Mm * ((m/m0)^(2/3)) * (1 - SRmin)
            res += f * kappa_diss * ssa * nm * molar_mass * pow(nm / nm0, 2.0 / 3.0) * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;

        }
        else // Omega > 1, precipitation kinetics
        {
            // Interfacial energy in J m-2 # surf_energy = 0.047
            const auto surf_energy = 0.047;

            // Molecular volume in m3 # molvol = 6.13E-29
            const auto molecular_volume = 6.13E-29;

            // u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 )
            const auto u = 16.0 * 3.14159 * std::pow(surf_energy, 3.0) * std::pow(molecular_volume, 2.0) / (3.0 * std::pow(1.38E-23 * T.val, 3.0));

            // Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009) # J0 = 1E20
            const auto j0 = 1E20;

            // Critical saturation threshold # SRcrit = exp( sqrt( u / LOG( J0 ) ) )
            const auto Omega_crit = std::exp(sqrt(u / std::log10(j0)));

            // Nucleation occurs only if saturation index above the threshold
            if(Omega >= Omega_crit)
            {
                // Nucleation rate in nuclie sec-1 # nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * Time
                const auto nuclie = j0 * exp(- u / std::pow(std::log10(Omega), 2.0));

                if(nuclie * dt >= 1) // Condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
                {
                    // Number of growth units in critical nuclie radius // nj = ( 2 * u ) / (LOG(SRmin)^3)
                    const auto nj = 2 * u / std::pow(std::log10(Omega), 3.0);

                    // Critical nuclie volume # vol = nj * molvol
                    const auto vol = nj * molecular_volume;

                    // Moles of nuclie formed # moles_nuc = ( -nuclie * vol ) / Mv
                    res_nuc += - nuclie * vol / molar_volume;
                }
            }

            // Neutral mechanism # knu = 1.8E-7 * exp((-66000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 1.8E-7 * exp(- 66.0 / R * (1.0 / T - 1.0 / 298.15));

            // Acidic mechanism # k1 = 1.9E-3 * exp((-67000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HCO3-") ^ 1.63)
            const auto kappa_acid = 1.9E-3 * exp(- 67.0 / R * (1.0 / T - 1.0 / 298.15)) * std::pow(activity_hco3, 1.63);

            // kpre = - (knu + k1)
            const auto kappa_pre = - (kappa_neu + kappa_acid);

            // moles_growth = -1 * S * m * Mm * kpre * (( ( SRmin^0.5 ) - 1 )^1.12) * Time
            res_growth += f * ssa * nm * molar_mass * kappa_pre * std::pow(std::pow(Omega, 0.5) - 1, 1.12);

            // Precipitation kinetic includes kinetic growth and nucleation
            res += res_growth + res_nuc;

            // TODO: implement upper bound in precipitation kinetics
            // 360 maxMol = MOL("Ca+2")
            // 370 IF (maxMol > MOL("HCO3-")) THEN maxMol = MOL("HCO3-")
            // 380 IF (maxMol < -moles) THEN moles = -maxMol

            // auto max_mol = mol_ca;
            // if(max_mol > mol_hco3) max_mol = mol_hco3;
            // if(max_mol < -res) res = - max_mol;
        }

        return res;
    };
    reaction_calcite.setRate(rate_func_calcite_olimse_dat);

    //Dolomite
    //CaMg(CO3)2 + 2.000H+ = 2.000HCO3- + 1.000Ca+2 + 1.000Mg+2
    //log_k     3.533
    //delta_h  -65.360    #kJ/mol        #95rob/hem
    //                                    -analytic -1.7923634E+3  -2.8963524E-1  9.9594493E+4  6.5114488E+2  -5.6008392E+6
    //#References = LogK/DGf: Internal calculation; DHf/DHr: 95rob/hem; S°: 95rob/hem; Cp: 95rob/hem; V°: 78hel/del,92ajoh;
    std::string eq_str_dolomite = "Dolomite + 2*H+ = 2*HCO3- + Ca++ + Mg++";
    MineralReaction min_reaction_dolomite = editor.addMineralReaction("Dolomite")
            .setEquation(eq_str_dolomite)
            .setSpecificSurfaceArea(0.13, "m2/g");
    Reaction reaction_dolomite = createReaction(min_reaction_dolomite, system);
    reaction_dolomite.setName("Dolomite reaction");
    ReactionRateFunction rate_func_dolomite_olimse_dat = [&min_reaction_dolomite, &reaction_dolomite, &system, &num_species, &R, &dt](const ChemicalProperties& properties) -> ChemicalScalar {

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0), res_growth(num_species, 0.0), res_nuc(num_species, 0.0);

        // Auxiliary variables with chemical scalars of value one
        ChemicalScalar f(num_species, 1.0);

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();
        const Pressure P = properties.pressure();
        const double P_atm = P.val * 1e-5 * 0.986923;

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_dolomite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_dolomite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0 = reaction_dolomite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_dolomite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        // Calculate TOT("Mg")
        // const Index i_mg = system.indexElementWithError("Mg");
        // const auto tot_mg = b(i_mg);
        // Calculate TOT("Ca")
        // const Index i_ca = system.indexElementWithError("Ca");
        // const auto tot_ca = b(i_ca);
        // Calculate TOT("C")
        // const Index i_c = system.indexElementWithError("c");
        // const auto tot_c = b(i_c);

        // Calculate activities ofr H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        //const Index i_co2g = system.indexSpeciesWithError("CO2(g)");
        //double activity_co2g = std::exp(lna(i_co2g));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 350
        // If (SRmin = 1) Then GoTo 350
        if((nm <= 0 && Omega < 1) || Omega == 1) // the is no way to dissolve further
            return res;

        // Average BET
        const auto ssa = min_reaction_dolomite.specificSurfaceArea();

        // Molar volume in m3/mol # Mv = 6.429E-05
        const auto molar_volume = 6.429E-05;

        if(Omega < 1) // Omega < 1, dissolution kinetics # kinetic data extracted from 04pal/kha
        {
            // Neutral mechanism # knu = 0.0000000295 * exp((-52200 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 0.0000000295 * exp(- 52.2 / R * (1.0 / T - 1.0 / 298.15));

            // Acidic mechanism # k1 = 0.000646 * exp((-36100 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.5)
            const auto kappa_acid = 0.000646 * exp(- 36.10 / R * (1.0 / T - 1.0 / 298.15)) * pow(activity_h, 0.5) ;

            // CO2(g) mechanism # k2 = 0.00000776 * exp((-34800 / 8.314) * ((1 / TK) - (1 / 298.15))) * (PR_P("CO2(g)") ^ 0.5)
            const auto kappa_co2g = 0.000646 * exp(- 34.80 / R * (1.0 / T - 1.0 / 298.15)) * pow(P_atm, 0.5) ;

            // total rate constant
            const auto kappa_diss = kappa_neu + kappa_acid + kappa_co2g;

            // Calculate the resulting mechanism function: rate = S * m * Mm * k * (1 - SRmin)
            res += f * kappa_diss * ssa * nm * molar_mass * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;

        }
        else // Omega > 1, precipitation kinetics
        {
            // Interfacial energy in J m-2 # surf_energy = 0.075
            const auto surf_energy = 0.075;

            // Molecular volume in m3 # molvol = 1.0674e-28
            const auto molecular_volume = 1.0674e-28;

            // u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 )
            const auto u = 16.0 * 3.14159 * std::pow(surf_energy, 3.0) * std::pow(molecular_volume, 2.0) / (3.0 * std::pow(1.38E-23 * T.val, 3.0));

            // Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009) # J0 = 1E20
            const auto j0 = 1E20;

            // Critical saturation threshold # SRcrit = exp( sqrt( u / LOG( J0 ) ) )
            const auto Omega_crit = std::exp(sqrt(u / std::log10(j0)));

            // Nucleation occurs only if saturation index above the threshold
            if(Omega >= Omega_crit)
            {
                // Nucleation rate in nuclie sec-1 # nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * time
                const auto nuclie = j0 * exp(- u / std::pow(std::log10(Omega), 2.0));

                if(nuclie * dt >= 1) // Condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
                {
                    // Number of growth units in critical nuclie radius // nj = ( 2 * u ) / (LOG(SRmin)^3)
                    const auto nj = 2 * u / std::pow(std::log10(Omega), 3.0);

                    // Critical nuclie volume # vol = nj * molvol
                    const auto vol = nj * molecular_volume;

                    // Moles of nuclie formed # moles_nuc = ( -nuclie * vol ) / Mv
                    res_nuc += - nuclie * vol / molar_volume;
                }
            }

            // Neutral mechanism # knu = 9.45E-15 * exp((-103100 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 9.45E-15 * exp(- 103.1 / R * (1.0 / T - 1.0 / 298.15));

            // kpre = - (knu + k1)
            const auto kappa_pre = - kappa_neu;

            // moles_growth = S * m * Mm * kpre * ABS(SRmin - 1) * Time
            res_growth += f * ssa * nm * molar_mass * kappa_pre * std::abs(Omega - 1);

            // Precipitation kinetic includes kinetic growth and nucleation
            res += res_growth + res_nuc;

            // TODO: implement upper bound in precipitation kinetics
            // 310 maxMol = TOT("Mg")
            // 320 IF (maxMol > TOT("Ca")) THEN maxMol = TOT("Ca")
            // 340 IF (maxMol > TOT("C")/2) THEN maxMol = TOT("C")/2
            // 350 IF (maxMol < -moles) THEN moles = -maxMol

            // auto max_mol = mol_mg;
            // if(max_mol > mol_ca) max_mol = mol_ca;
            // if(max_mol > mol_c / 2) max_mol = mol_c / 2;
            // if(max_mol < -res) res = - max_mol;
        }

        return res;

    };
    reaction_dolomite.setRate(rate_func_dolomite_olimse_dat);

    // Create the ReactionSystem instances
    ReactionSystem reactions(system, {reaction_calcite, reaction_dolomite});

    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Calcite", "Dolomite"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("O2",    1.0, "umol");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("SiO2",  10,  "mol");
    problem_ic.add("MgCl2", 1e-10, "mol");

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O",   1.00, "kg");
    problem_bc.add("O2",    1.0, "umol");
    problem_bc.add("NaCl",  0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2",   0.75, "mol");

    //equilibrium_options.optimum.output.active = true;
    //ChemicalState state_ic = equilibrate(problem_ic, equilibrium_options);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    //std::cout << "state_ic = " << state_ic << std:: endl;
    //std::cout << "state_bc = " <<state_bc << std:: endl;

    //state_ic.setSpeciesAmount("Dolomite", 10, "mcmol");

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.81, "m3");   // 0.81 = 0.90 * 0.9 (0.9 of rock is due to 10% porosity, 0.90 is 90% quartz of the rock)
    state_ic.scalePhaseVolume("Calcite", 0.054, "m3");  // 0.054 = 0.06 * 0.9 (0.9 of rock due to 10% porosity, 0.06 is 6% calcite of the rock)
    state_ic.scalePhaseVolume("Dolomite", 0.036, "m3");  // 0.036 = 0.04 * 0.9 (0.9 of rock due to 10% porosity, 0.04 is 4% dolomite of the rock)

    // Set initial value of minerals reacting
    reaction_calcite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_dolomite.setInitialAmounts(state_ic.speciesAmounts());

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    reactive_transport_options.use_smart_kinetic_solver = params.use_smart_kinetics_solver;

    reactive_transport_options.equilibrium = equilibrium_options;
    reactive_transport_options.smart_equilibrium = smart_equilibrium_options;
    reactive_transport_options.kinetics = kinetic_options;
    reactive_transport_options.smart_kinetics = smart_kinetic_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(reactions, partition);
    rtsolver.setOptions(reactive_transport_options);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(params.v);
    rtsolver.setDiffusionCoeff(params.D);
    rtsolver.setBoundaryState(state_bc);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Ca++)");
    output.add("speciesMolality(Mg++)");
    output.add("speciesMolality(HCO3-)");
    output.add("speciesMolality(CO2(aq))");
    output.add("phaseVolume(Calcite)");
    output.add("phaseVolume(Dolomite)");
    output.add("speciesMolality(CO3--)");
    output.add("speciesMolality(CaCl+)");
    output.add("speciesMolality(Ca(HCO3)+)");
    output.add("speciesMolality(MgCl+)");
    output.add("speciesMolality(Mg(HCO3)+)");
    output.add("speciesMolality(OH-)");
    output.add("elementmolality(C)");
    output.add("elementmolality(Ca)");
    output.add("elementmolality(Cl)");
    output.add("elementmolality(H)");
    output.add("elementmolality(Mg)");
    output.add("elementmolality(Na)");
    output.add("elementmolality(O)");
    output.add("elementmolality(Si)");
    output.add("elementmolality(Z)");
    output.add("speciesMolality(MgCO3(aq))");
    output.add("speciesMolality(MgOH+)");
    output.add("speciesAmount(Calcite)");
    output.add("speciesAmount(Dolomite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(TRANSPORT)

    // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print some progress
        //std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        ReactiveTransportResult result = rtsolver.result();

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        results.time_reactive_transport_smart_kin_smart_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        results.time_reactive_transport_conv_kin_conv_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        results.time_reactive_transport_conv_kin_smart_eq = toc(TRANSPORT);
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        results.time_reactive_transport_smart_kin_conv_eq = toc(TRANSPORT);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-smart-kin-smart-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-conventional-kin-conventional-eq.json") << analysis;
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-conventional-kin-smart-eq.json") << analysis;
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver)
        JsonOutput(folder + "/" + "analysis-smart-kin-conventional-eq.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) {
        results.smart_kinetic_timing = analysis.smart_kinetics.timing;
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_kin_smart_eq_acceptance_rate = analysis.smart_kinetics.smart_kinetics_estimate_acceptance_rate;
        results.smart_kin_smart_eq_equilibrium_acceptance_rate
                = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
        results.outputStatisticsSmartKineticsSmartEquilibrium(params);
    }
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) {
        results.kinetic_timing = analysis.kinetics.timing;
        results.equilibrium_timing = analysis.equilibrium.timing;
        results.outputStatisticsConventionalKineticsConventionalEquilibrium(params);
    }
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.kinetic_timing = analysis.kinetics.timing;
        results.conv_kin_smart_eq_equilibrium_acceptance_rate
                = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
        results.outputStatisticsConventionalKineticsSmartEquilibrium(params);
    }
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) {
        results.smart_kinetic_timing = analysis.smart_kinetics.timing;
        results.equilibrium_timing = analysis.equilibrium.timing;
        results.smart_kin_conv_eq_acceptance_rate = analysis.smart_kinetics.smart_kinetics_estimate_acceptance_rate;
        results.outputStatisticsSmartKineticsConventionalEquilibrium(params);
    }
}

