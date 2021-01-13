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
    int minute = 60;
    int hour = 60 * minute;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 25.0; // the x-coordinates of the right boundaries
    params.ncells = 243; // the number of cells in the spacial discretization
    params.nsteps = 1000; // the total number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 1*hour; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;             // the diffusion coefficient (in units of m2/s)
    params.v = 0.8e-5;          // the Darcy velocity (in units of m/s)
    params.T = 60.0;            // the temperature (in units of degC)
    params.P = 200 * 1.01325;   // the pressure (in units of bar)

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // *********************************************************************************************************** //
    // Clustering approach
    // *********************************************************************************************************** //

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

//    // *********************************************************************************************************** //
//    // Priority queue approach
//    // *********************************************************************************************************** //
//
//    // ----------------------------------------------------------- //
//    // Smart equilibrium parameters
//    // ----------------------------------------------------------- //
//
//    // Run smart algorithm with clustering
//    params.method = SmartEquilibriumStrategy::PriorityQueue;
//    params.smart_equilibrium_reltol = 1e-2;
//
//    // ----------------------------------------------------------- //
//    // Smart kinetic parameters
//    // ----------------------------------------------------------- //
//
//    // Select priority queue approach
//    params.smart_kin_method = SmartKineticStrategy::PriorityQueue;
//    params.smart_kinetic_tol = 1e-2;
//    params.smart_kinetic_abstol = 1e-5;
//    params.smart_kinetic_reltol = 1e-1;

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
    /// Execute reactive transport with different solvers
    //params.use_smart_kinetics_solver = true; params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

    /// **************************************************************************************************************///
    /// CONVENTIONAL kinetics & CONVENTIONAL equilibrium
    /// **************************************************************************************************************///
    params.use_smart_kinetics_solver = false; params.use_smart_equilibrium_solver = false;
    params.outputConsoleKineticMethod();
    runReactiveTransport(params, results);

    /// **************************************************************************************************************///
    /// SMART kinetics & SMART equilibrium
    /// **************************************************************************************************************///
    params.use_smart_kinetics_solver = true;  params.use_smart_equilibrium_solver = true;
    params.outputConsoleKineticMethod();
    runReactiveTransport(params, results);

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
    auto folder = params.makeResultsFolderKinetics("results-kinetics-scaling");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    //equilibrium_options.hessian = params.hessian;

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
    StringList selected_elements = "H Cl S O Ba Ca Sr Na K Mg C Si";

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == ActivityModel::HKF){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == ActivityModel::Pitzer){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ActivityModel::DebyeHuckel){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Barite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Create reaction for the Barite mineral
    std::string eq_str_barite = "Barite = SO4-- + Ba++";
    MineralReaction min_reaction_barite = editor.addMineralReaction("Barite")
            .setEquation(eq_str_barite)
            .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol")
            .setSpecificSurfaceArea(0.006, "m2/g");
    Reaction reaction_barite = createReaction(min_reaction_barite, system);
    reaction_barite.setName("Barite kinetics");

    // Ionic strength function
    const auto I = ChemicalProperty::ionicStrength(system);
    // Barite kinetic rate function
    ReactionRateFunction rate_func_barite_shell = [&min_reaction_barite, &reaction_barite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_barite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_barite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0 = reaction_barite.initialAmounts();

//        // Amount of elements
//        const auto b = system.elementAmounts(n.val);
//
//        // Calculate TOT("Ba")
//        const Index i_ba = system.indexElementWithError("Ba");
//        const auto b_ba = b(i_ba);
//        const auto tot_ba = b_ba;
//
//        // Calculate TOT("S(6)") = TOT("SO4")
//        const Index i_s = system.indexElementWithError("S");
//        const auto b_s = b(i_s);
//        const Index i_o = system.indexElementWithError("O");
//        const auto b_o = b(i_o);
//        const auto tot_so4 = b_s + 4 * b_o;

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_barite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        // Prevent negative mole numbers here for the solution of the ODEs
        nm = std::max(nm, 0.0);
        nm0 = std::max(nm0, 0.0);

        // Get H+ activity and ionic strength
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        double ionic_strength = I(properties).val;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 240
        // If (SRmin = 1) Then GoTo 240
        if((nm <= 0 && Omega < 1) || (Omega == 1)) // the is no way to precipitate further
            return res;
        // S = 0.006 # average BET from 16zhe/did ; suggested value in m2/g
        //const auto ssa = 0.006 * 1e3;
        const auto ssa = 0.006;

        // If (SRmin > 1) Then GoTo 130
        if(Omega > 1) // precipitation kinetics
        {
            // If (m <= 1e-5) then GoTo 170
            if(nm > 1e-5)
            {
                // knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
                // kpre = (-1) * knu
                const auto kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15))
                                       * pow(10, 0.6 * pow(ionic_strength, 0.5));

                // rate = S * m * Mm * kpre * (SRmin - 1)
                res += f * ssa * nm * molar_mass * kappa_pre * (Omega - 1);
            }
            else
                // Set nucleation rate
                res = -1e-12 * f;

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_ba;
//            if(max_mol > tot_so4) max_mol = tot_so4;
//            if(max_mol < -res) return -max_mol * f;


        }
        else // dissolution kinetics
        {
            // k1 = 2.75E-08 * exp((-25000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.03) * (10 ^ (0.6 * MU^0.5))
            const auto kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15))
                               * std::pow(activity_h, 0.03)
                               * pow(10, 0.6 * pow(ionic_strength, 0.5));

            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - (SRmin ^ 0.2))
            res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - pow(Omega, 0.2));

//            // Do not dissolve more than what is available
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res +=  nm.val;
        }
        return res;

    };
    reaction_barite.setRate(rate_func_barite_shell);

    // Step **: Create the ReactionSystem instances
    //ReactionSystem reactions(editor);
    ReactionSystem reactions(system, {reaction_barite});

    // Step **: Create the ReactionSystem instances
    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Barite"});

    // ************************************************************************************************************** //
    // Initial condition (IC)
    // ************************************************************************************************************** //

    double water_kg = 1.0;
    // Define the initial condition *with formation water)
    EquilibriumInverseProblem problem_ic_fw(system);
    problem_ic_fw.setTemperature(params.T, "celsius");
    problem_ic_fw.setPressure(params.P, "atm");
    problem_ic_fw.add("H2O", water_kg, "kg");
    problem_ic_fw.add("SO4", 10 * water_kg, "ug");
    problem_ic_fw.add("Ca", 995 * water_kg, "mg");
    problem_ic_fw.add("Ba", 995 * water_kg, "mg");
    problem_ic_fw.add("Sr", 105 * water_kg, "mg");
    problem_ic_fw.add("Na", 27250 * water_kg, "mg");
    problem_ic_fw.add("K", 1730 * water_kg, "mg");
    problem_ic_fw.add("Mg", 110 * water_kg, "mg");
    problem_ic_fw.add("Cl", 45150 * water_kg, "mg");
    problem_ic_fw.add("HCO3", 1980 * water_kg, "mg");
    problem_ic_fw.pH(7.0, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_ic = equilibrate(problem_ic_fw);

    //state_ic.setSpeciesAmount("Barite", 1.0, "umol");
    state_ic.setSpeciesAmount("Barite", 10, "mcmol");

    // Scale the percentage of the aqueous phase and the whole volume
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // Set initial value of Barite
    reaction_barite.setInitialAmounts(state_ic.speciesAmounts());

    // Define function to evaluate ph of the chemical system
    auto evaluate_pH = ChemicalProperty::pH(system);

    // Fetch the pH of the initial sate
    ChemicalProperties props = state_ic.properties();
    ChemicalScalar pH_ic = evaluate_pH(props);
    std::cout << "ph(FW) = " << pH_ic.val << std:: endl;
    //std::cout << "state_ic b: " << state_ic << std::endl;
    //getchar();

    // ************************************************************************************************************** //
    // Boundary condition (BC)
    // ************************************************************************************************************** //

    // Define the first boundary condition (with completion brine)
    EquilibriumProblem problem_bc_cb(system);
    problem_bc_cb.setTemperature(params.T, "celsius");
    problem_bc_cb.setPressure(params.P, "atm");
    problem_bc_cb.add("H2O", water_kg, "kg");
    problem_bc_cb.add("NaCl", 0.7, "mol");

    // Equilibrate the initial condition
    ChemicalState state_bc_cb = equilibrate(problem_bc_cb);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_cb.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_cb.properties();
    ChemicalScalar pH_bc = evaluate_pH(props);
    std::cout << "ph(CB) = " << pH_bc.val << std:: endl;

    // Define the first boundary condition (with seawater)
    EquilibriumInverseProblem problem_bc_sw(system);
    problem_bc_sw.setTemperature(params.T, "celsius");
    problem_bc_sw.setPressure(params.P, "atm");
    problem_bc_sw.add("H2O", water_kg, "kg");
    problem_bc_sw.add("SO4--", 2710 * water_kg, "mg");
    problem_bc_sw.add("Ca++", 411 * water_kg, "mg");
    problem_bc_sw.add("Ba++", 0.01 * water_kg, "mg");
    problem_bc_sw.add("Sr++", 8 * water_kg, "mg");
    problem_bc_sw.add("Na+", 10760 * water_kg, "mg");
    problem_bc_sw.add("K+", 399 * water_kg, "mg");
    problem_bc_sw.add("Mg++", 1290 * water_kg, "mg");
    problem_bc_sw.add("Cl-", 19350 * water_kg, "mg");
    problem_bc_sw.add("HCO3-", 142 * water_kg, "mg");
    problem_bc_sw.pH(8.1, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_bc_sw = equilibrate(problem_bc_sw);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_sw.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_sw.properties();
    ChemicalScalar pH_sw = evaluate_pH(props);
    std::cout << "ph(SW) = " << pH_sw.val << std:: endl;

    //std::cout << "state_bc_sw b: " << tr(state_bc_sw.elementAmounts()) << std::endl;
    //getchar();

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
    rtsolver.setBoundaryState(state_bc_sw);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    if(params.output_results)
    {
        ChemicalOutput output(rtsolver.output());
        output.add("pH");
        output.add("speciesAmount(H+)");
        output.add("speciesAmount(Cl-)");
        output.add("speciesAmount(SO4--)");
        output.add("speciesAmount(Ba++)");
        output.add("speciesAmount(Ca++)");
        output.add("speciesAmount(Sr++)");
        output.add("speciesAmount(Na+)");
        output.add("speciesAmount(Barite)");
        output.add("elementAmount(Ba)");
        output.add("elementAmount(C)");
        output.add("elementAmount(Ca)");
        output.add("elementAmount(Cl)");
        output.add("elementAmount(H)");
        output.add("elementAmount(K)");
        output.add("elementAmount(Mg)");
        output.add("elementAmount(Na)");
        output.add("elementAmount(O)");
        output.add("elementAmount(S)");
        output.add("elementAmount(Si)");
        output.add("elementAmount(Sr)");
        output.add("elementAmount(Z)");
        output.add("phaseAmount(Barite)");
        output.add("phaseMass(Barite)");
        output.add("phaseVolume(Barite)");

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
        //if (!(step % 100))
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

