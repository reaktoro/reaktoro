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
    Time start = time();

    // Step 1: Initialise auxiliary time-related constants
    // int second = 1;
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0;                                        // the x-coordinates of the left boundaries
    params.xr = 100.0;                                      // the x-coordinates of the right boundaries
    params.ncells = 100;                                    // the number of cells in the spacial discretization
    params.nsteps = 1000;                                   // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells;    // the time step (in units of s)
    params.dt = 0.1*day;                                    // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 3e-5;    // the Darcy velocity (in units of m/s)
    params.T = 25.0;    // the temperature (in units of degC)
    params.P = 1.01325; // the pressure (in units of bar)


    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

//    // *********************************************************************************************************** //
//    // Clustering approach
//    // *********************************************************************************************************** //

//    // ----------------------------------------------------------- //
//    // Smart equilibrium parameters
//    // ----------------------------------------------------------- //
//
//    // Run smart algorithm with clustering
//    params.method = SmartEquilibriumStrategy::Clustering;
//    params.smart_equilibrium_reltol = 1e-3;
//
//    // ----------------------------------------------------------- //
//    // Smart kinetic parameters
//    // ----------------------------------------------------------- //
//
//    // Select clustering approach
//    params.smart_kin_method = SmartKineticStrategy::Clustering;
//    params.smart_kinetic_tol = 1e-3; // 5e-3; // 1e-4;
//    params.smart_kinetic_reltol = 1e-1;
//    params.smart_kinetic_abstol = 1e-4;

    // *********************************************************************************************************** //
    // Priority queue approach
    // *********************************************************************************************************** //

    // ----------------------------------------------------------- //
    // Smart equilibrium parameters
    // ----------------------------------------------------------- //

    // Run smart algorithm with clustering
    params.method = SmartEquilibriumStrategy::PriorityQueue;
    params.smart_equilibrium_reltol = 1e-2;

    // ----------------------------------------------------------- //
    // Smart kinetic parameters
    // ----------------------------------------------------------- //

    // Select priority queue approach
    params.smart_kin_method = SmartKineticStrategy::PriorityQueue;
    params.smart_kinetic_tol = 1e-2;
    params.smart_kinetic_abstol = 1e-5;
    params.smart_kinetic_reltol = 1e-1;

//    // *********************************************************************************************************** //
//    // Nearest neighbour approach
//    // *********************************************************************************************************** //

//    // ----------------------------------------------------------- //
//    // Smart equilibrium parameters
//    // ----------------------------------------------------------- //

//    // Run smart algorithm with clustering
//    params.method = SmartEquilibriumStrategy::NearestNeighbour;
//    params.smart_equilibrium_reltol = 1e-1;
//
//    // ----------------------------------------------------------- //
//    // Smart kinetic parameters
//    // ----------------------------------------------------------- //
//
//    // Select nearest neighbour approach
//    params.smart_kin_method = SmartKineticStrategy::NearestNeighbour;
//    params.smart_kinetic_abstol = 1e-4;
//    params.smart_kinetic_reltol = 1e-1;

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

    /// **************************************************************************************************************///
    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
    /// **************************************************************************************************************///
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
    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportKineticsResults& results) -> void
{
    // Step **: Create the results folder
    auto folder = params.makeResultsFolderKinetics("scavenging-hematite");

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
    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);

    StringList selected_species = {"Ca(HCO3)+", "CO3--", "CaCO3(aq)", "Ca++", "CaSO4(aq)", "CaOH+", "Cl-",
                                   "FeCl++", "FeCl2(aq)", "FeCl+", "Fe++", "FeOH+",  "FeOH++", "Fe+++",
                                   "H2(aq)", "HSO4-", "H2S(aq)", "HS-", "H2O(l)",  "H+", "OH-", "HCO3-",
                                   "K+", "KSO4-", "Mg++", "MgSO4(aq)", "MgCO3(aq)", "MgOH+", "Mg(HCO3)+",
                                   "Na+", "NaSO4-", "O2(aq)", "S5--", "S4--", "S3--", "S2--", "SO4--"};
    StringList selected_elements = "C Ca Cl Fe H K Mg Na O S";

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
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Pyrite");
    editor.addMineralPhase("Hematite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //std::cout << "system = \n" << system << std:: endl;

    // Hematite: Fe2O3
    std::string eq_str_hematite = "Hematite + 6*H+ = 3*H2O(l) + 2*Fe+++";
    MineralReaction min_reaction_hematite = editor.addMineralReaction("Hematite")
            .setEquation(eq_str_hematite)
            .setSpecificSurfaceArea(1.0, "m2/g");
    Reaction reaction_hematite = createReaction(min_reaction_hematite, system);
    reaction_hematite.setName("Hematite reaction");
    ReactionRateFunction rate_func_hematite_shell = [&min_reaction_hematite, &reaction_hematite, &system](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_hematite.equation(), system);

        // Auxiliary variables
        ChemicalScalar f(num_species, 1.0);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_hematite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_hematite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition().val;
        const auto n0 = reaction_hematite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_hematite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral];
        auto nm0 = n0[imineral];


        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_hs = system.indexSpeciesWithError("HS-");
        double activity_hs = std::exp(lna(i_hs));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        /*
         * Hematite
        # kinetic data extracted from 04pal/kha 04pou/kro
        # warning dissolution only
        # Confidence level: 4
        -start
        1 SRmin = SR("Hematite")
        10 moles = 0
        20 If (m <= 0) and (SRmin < 1) Then GoTo 250
        30 S = PARM(1) # Default value
        40 Mm = 159.7 # molar mass in g/mol
        50 If (SRmin > 1) Then GoTo 250
        ########## start dissolution bloc ##########
        60 knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
        70 k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
        80 k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
        90 k = knu + k1 + k2
        100 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
        110 moles = rate * Time
        120 REM Do not dissolve more than what is available
        130 IF (moles > M) THEN moles = M
        ########## end dissolution bloc ##########
        250 Save moles
        -end
         */
        // If (m <= 0) and (SRmin < 1) Then GoTo 250
        if(nm <= 0 && Omega < 1) // the is no way to precipitate further
            res = ChemicalScalar(num_species, 0.0);
        // S = PARM(1) # Default value
        const auto ssa = min_reaction_hematite.specificSurfaceArea();

        if(Omega < 1) // dissolution kinetics
        {
            /*
             * ########## start dissolution bloc ##########
            60 knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
            70 k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
            80 k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
            90 k = knu + k1 + k2
            100 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            110 moles = rate * Time
            120 REM Do not dissolve more than what is available
            130 IF (moles > M) THEN moles = M
            ########## end dissolution bloc ##########
             */

            // knu = 2.51E-15 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_1 = 2.51e-15 * exp(- 66.2 / R * (1.0/T - 1.0/298.15));

            // k1 = 0.000000000407 * exp((-66200 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 1.0)
            const auto kappa_2 = 0.000000000407 * exp( -66.2 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 1.0);

            // Sulfide catalyzer (sulfide promotion)
            // k2 = 3.50E-9 * exp((-40000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("HS-") ^ 0.5)
            const auto kappa_3 = 3.50e-9 * exp(- 40 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_hs, 0.5);

            const auto kappa = kappa_1 + kappa_2 + kappa_3;
            //std::cout << "kappa = " << kappa << std:: endl;

            // Calculate the resulting mechanism function
            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm * molar_mass * pow(nm/nm0, 2/3) * kappa * (1 - Omega);

            // Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (res > nm.val) res += nm.val;

        }

        return res;

    };
    reaction_hematite.setRate(rate_func_hematite_shell);

    // Step **: Create the ReactionSystem instances
    ReactionSystem reactions(system, {reaction_hematite});

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Hematite"});

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumInverseProblem problem_ic(partition);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 58.0, "kg");
    problem_ic.add("Cl-", 1122.3e-3, "kg");
    problem_ic.add("Na+", 624.08e-3, "kg");
    problem_ic.add("SO4--", 157.18e-3, "kg");
    problem_ic.add("Mg++", 74.820e-3, "kg");
    problem_ic.add("Ca++", 23.838e-3, "kg");
    problem_ic.add("K+", 23.142e-3, "kg");
    problem_ic.add("HCO3-", 8.236e-3, "kg");
    problem_ic.add("O2(aq)", 58e-12, "kg");
    problem_ic.pH(8.951);
    problem_ic.pE(8.676);

    ChemicalState state_ic = equilibrate(problem_ic);

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumInverseProblem problem_bc(partition);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 58.0, "kg");
    problem_bc.add("Cl-", 1122.3e-3, "kg");
    problem_bc.add("Na+", 624.08e-3, "kg");
    problem_bc.add("SO4--", 157.18e-3, "kg");
    problem_bc.add("Mg++", 74.820e-3, "kg");
    problem_bc.add("Ca++", 23.838e-3, "kg");
    problem_bc.add("K+", 23.142e-3, "kg");
    problem_bc.add("HCO3-", 8.236e-3, "kg");
    problem_bc.add("O2(aq)", 58e-12, "kg");
    problem_bc.add("HS-", 0.0196504, "mol");
    problem_bc.add("H2S(aq)", 0.167794, "mol");
    problem_bc.pH(5.726);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Set initial Hematite amount
    state_ic.setSpeciesAmount("Hematite", 0.5, "mol"); // MM(Hematite) = 115.86 g / mol

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // Set initial value of Hemitite
    reaction_hematite.setInitialAmounts(state_ic.speciesAmounts());

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
    if(params.output_results)
    {
        ChemicalOutput output(rtsolver.output());
        output.add("pH");
        output.add("speciesAmount(H+)");
        output.add("speciesAmount(HS-)");
        output.add("speciesAmount(S2--)");
        output.add("speciesAmount(CO3--)");
        output.add("speciesAmount(HSO4-)");
        output.add("speciesAmount(H2S(aq))");
        output.add("speciesAmount(Fe++)");
        output.add("speciesAmount(Fe+++)");
        output.add("speciesAmount(Siderite)");
        output.add("speciesAmount(Pyrite)");
        output.add("speciesAmount(Hematite)");
        output.add("phaseAmount(Siderite)");
        output.add("phaseAmount(Pyrite)");
        output.add("phaseAmount(Hematite)");
        output.add("phaseMass(Siderite)");
        output.add("phaseMass(Pyrite)");
        output.add("phaseMass(Hematite)");
        output.filename(folder + "/" + "test.txt");
    }

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(TRANSPORT);

    // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print some progress
        //if (!(step % 1))
        //std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.time_reactive_transport_smart_kin_smart_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.time_reactive_transport_conv_kin_conv_eq = toc(TRANSPORT);
    if(!params.use_smart_kinetics_solver && params.use_smart_equilibrium_solver) results.time_reactive_transport_conv_kin_smart_eq = toc(TRANSPORT);
    if(params.use_smart_kinetics_solver && !params.use_smart_equilibrium_solver) results.time_reactive_transport_smart_kin_conv_eq = toc(TRANSPORT);

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