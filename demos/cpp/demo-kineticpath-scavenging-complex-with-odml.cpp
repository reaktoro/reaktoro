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

using namespace Reaktoro;

// Reactive transport test includes
#include <demos/cpp/DemosUtils.h>

auto runKinetics(KineticPathParams &) -> void;

int main()
{
    int minute = 60;
    int hour = 60 * minute;

    KineticPathParams params = {};

    params.t0 = 0;

    params.dt = hour;
    params.n = 1400;
    //params.n = 10;

    params.tfinal = params.n * hour;

    params.T = 25;
    params.P = 1.01325; // 1 atm = 1.01325 bar

    params.activity_model = ActivityModel::DebyeHuckel;

    // ----------------------------------------------------------- //
    // Smart equilibrium parameters
    // ----------------------------------------------------------- //

    params.use_smart_equilibrium_solver = false;

    params.smart_eq_method = SmartEquilibriumStrategy::Clustering;
    params.smart_equilibrium_reltol = 1e-1;

    //params.method = SmartEquilibriumStrategy::PriorityQueue;
    //params.smart_equilibrium_reltol = 1e-1;

    //params.method = SmartEquilibriumStrategy::NearestNeighbour;
    //params.smart_equilibrium_reltol = 1e-2;

    // ----------------------------------------------------------- //
    // Smart kinetic parameters
    // ----------------------------------------------------------- //

    // Select clustering approach
    params.smart_kin_method = SmartKineticStrategy::Clustering;
    //params.smart_kinetic_tol = 1e-2;
    //params.smart_kinetic_tol = 1e-3;
    params.smart_kinetic_tol = 1e-4;
    //params.smart_kinetic_tol = 1e-5;
    //params.smart_kinetic_tol = 1e-6;
    params.smart_kinetic_abstol = 1e-4;
    params.smart_kinetic_reltol = 1e-1;

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

    params.use_smart_kinetic_solver = true; runKinetics(params);
    //params.use_smart_kinetic_solver = false; runKinetics(params);
}

auto runKinetics(KineticPathParams & params) -> void
{
    // Create results folder
    auto filename = params.makeResultsFile("scavenging-complex-no-nuc");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = GibbsHessian::Exact;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.learning.hessian = equilibrium_options.hessian;
    smart_equilibrium_options.method = params.smart_eq_method;

    // Step **: Define chemical kinetic solver options
    KineticOptions kinetic_options;
    kinetic_options.equilibrium = equilibrium_options;
    kinetic_options.smart_equilibrium = smart_equilibrium_options;
    kinetic_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;

    // Step **: Define smart chemical kinetic solver options
    SmartKineticOptions smart_kinetic_options;
    smart_kinetic_options.tol = params.smart_kinetic_tol;
    smart_kinetic_options.reltol = params.smart_kinetic_reltol;
    smart_kinetic_options.abstol = params.smart_kinetic_abstol;
    smart_kinetic_options.learning = kinetic_options;
    smart_kinetic_options.learning.equilibrium = equilibrium_options;
    smart_kinetic_options.smart_equilibrium = smart_equilibrium_options;
    smart_kinetic_options.use_smart_equilibrium_solver = kinetic_options.use_smart_equilibrium_solver;
    smart_kinetic_options.method = params.smart_kin_method;

    KineticPathOptions kineticpath_options;
    kineticpath_options.use_smart_kinetic_solver = params.use_smart_kinetic_solver;
    kineticpath_options.use_smart_equilibrium_solver = kinetic_options.use_smart_equilibrium_solver;
    kineticpath_options.equilibrium = equilibrium_options;
    kineticpath_options.smart_equilibrium = smart_equilibrium_options;
    kineticpath_options.kinetics = kinetic_options;
    kineticpath_options.smart_kinetics = smart_kinetic_options;

    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);
    StringList selected_elements = "Al C Ca Cl Fe H K Mg Na O S Si";

    if(params.activity_model == ActivityModel::Pitzer)
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    else if(params.activity_model == ActivityModel::DebyeHuckel)
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);

    editor.addMineralPhase("Pyrrhotite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Daphnite,14A"); // Chlorite
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Kaolinite");
    editor.addMineralPhase("Quartz");

    // Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

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
        if((nm <= 0 && Omega < 1) || Omega == 1) // the is no way to precipitate further
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

//            // Nucleation occurs only if saturation index above the threshold
//            if(Omega >= Omega_crit)
//            {
//                // Nucleation rate in nuclie sec-1 # nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * Time
//                const auto nuclie = j0 * exp(- u / std::pow(std::log10(Omega), 2.0));
//
//                if(nuclie * dt >= 1) // Condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
//                {
//                    // Number of growth units in critical nuclie radius // nj = ( 2 * u ) / (LOG(SRmin)^3)
//                    const auto nj = 2 * u / std::pow(std::log10(Omega), 3.0);
//
//                    // Critical nuclie volume # vol = nj * molvol
//                    const auto vol = nj * molecular_volume;
//
//                    // Moles of nuclie formed # moles_nuc = ( -nuclie * vol ) / Mv
//                    res_nuc += - nuclie * vol / molar_volume;
//                }
//            }

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

    //Chamosite(Daphnite)
    //Daphnite,14A in supcrt07 database
    //Fe5Al(AlSi3)O10(OH)8 + 16.000H+ = 2.000Al+3 + 5.000Fe+2 + 3.000H4SiO4 + 6.000H2O
    //     log_k    47.579
    //     delta_h -504.518    #kJ/mol        #01vid/par
    //     -analytic -2.6210061E+3  -4.2497094E-1  1.5576281E+5  9.4858884E+2  -6.610337E+6
    //     #References = LogK/DGf: Internal calculation; DHf/DHr: 01vid/par; S°: 01vid/par; Cp: 05vid/par; V°: 05vid/par;
    //Fe5Al(AlSi3)O10(OH)8 + 16.000H+ = 2.000Al+3 + 5.000Fe+2 + 3.000 (SiO2 + 2H2O) + 6.000H2O
    std::string eq_str_daphnite = "Daphnite,14A + 16*H+ = 2*Al+++ + 5*Fe++ + 3*SiO2(aq) + 12*H2O(l)";
    MineralReaction min_reaction_daphnite = editor.addMineralReaction("Daphnite,14A")
            .setEquation(eq_str_daphnite)
            .setSpecificSurfaceArea(0.0027, "m2/g"); //  S = 0.0027, 0.2% BET (03bra/bos) suggested value in m2/g
    Reaction reaction_daphnite = createReaction(min_reaction_daphnite, system);
    reaction_daphnite.setName("Daphnite reaction");
    ReactionRateFunction rate_func_daphnite_olimse_dat = [&min_reaction_daphnite, &reaction_daphnite, &system, &num_species, &R](const ChemicalProperties& properties) -> ChemicalScalar {
        
        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0);

        // Auxiliary variables with chemical scalars of value one
        ChemicalScalar f(num_species, 1.0);
        
        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_daphnite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_daphnite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0 = reaction_daphnite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_daphnite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        // Calculate activities for H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_oh = system.indexSpeciesWithError("OH-");
        double activity_oh = std::exp(lna(i_oh));

        // The molar mass of the mineral (in units of kg/mol) # Mm = 713.5 # molar mass in g/mol
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 250
        if(nm <= 0 && Omega < 1) // the is no way to precipitate further
            res = ChemicalScalar(num_species, 0.0);

        // 0.2% BET (03bra/bos)
        const auto ssa = min_reaction_daphnite.specificSurfaceArea();

        if(Omega < 1) // Omega < 1, dissolution kinetics
        {
            // Neutral mechanism: knu = 6.4E-17 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 6.4E-17 * exp(- 16.0 / R * (1.0/T - 1.0/298.15));

            // Acidic mechanism: k1 = 8.2E-09 * exp((-17000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H+") ^ 0.28)
            const auto kappa_acid = 8.2E-09 * exp(- 17.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.28);

            // Mechanism with OH- catalyzer (hydroxide promotion): k2 = 6.9E-09 * exp((-16000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.34)
            const auto kappa_oh = 6.9E-09 * exp(- 16.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_oh, 0.34);

            const auto kappa_diss = kappa_neu + kappa_acid + kappa_oh;

            // Calculate the resulting mechanism function: rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin)
            res += f * ssa * nm * molar_mass * pow(nm / nm0, 2.0 / 3.0) * kappa_diss * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;
        }

        return res;
    };
    reaction_daphnite.setRate(rate_func_daphnite_olimse_dat);

    //Siderite
    //FeCO3 + 1.000H+ = 1.000HCO3- + 1.000Fe+2
    //     log_k    -0.273
    //     delta_h  -27.862    #kJ/mol        #Internal calculation
    //     -analytic -9.0291123E+2  -1.4586221E-1  4.9931005E+4  3.2756219E+2  -2.8333834E+6
    //     #References = LogK/DGf: 04chi; DHf/DHr: Internal calculation; S°: 04chi; Cp: 04chi; V°: 78hel/del,85hel;
    std::string eq_str_siderite = "Siderite + H+ = HCO3- + Fe++";
    MineralReaction min_reaction_siderite = editor.addMineralReaction("Siderite")
            .setEquation(eq_str_siderite)
            .setSpecificSurfaceArea(0.13, "m2/g"); // S = 0.13, suggested value in m2/g
    Reaction reaction_siderite = createReaction(min_reaction_siderite, system);
    reaction_siderite.setName("Siderite reaction");
    ReactionRateFunction rate_func_siderite_olimse_dat = [&min_reaction_siderite, &reaction_siderite, &system, &num_species, &R, &dt](const ChemicalProperties& properties) -> ChemicalScalar {

        // The mineral reaction rate using specified surface area
        ChemicalScalar res(num_species, 0.0), res_growth(num_species, 0.0), res_nuc(num_species, 0.0);

        // Auxiliary variables with chemical scalars of value one
        ChemicalScalar f(num_species, 1.0);
        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_siderite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_siderite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        const auto n0 = reaction_siderite.initialAmounts();

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_siderite.mineral());

        // The number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        // Calculate TOT("Fe")
        // const Index i_fe = system.indexElementWithError("Fe");
        // const auto tot_fe = b(i_fe);
        // Calculate TOT("C")
        // const Index i_c = system.indexElementWithError("c");
        // const auto tot_c = b(i_c);

        // Calculate activities for H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 350
        // If (SRmin = 1) Then GoTo 350
        if((nm <= 0 && Omega < 1) || Omega == 1) // the is no way to precipitate further
            return res;

        // Average BET
        const auto ssa = min_reaction_siderite.specificSurfaceArea();

        // Molar mass, Mv = 2.049E-05, in g/mol
        const auto molar_volume = 2.049E-05;

        if(Omega < 1) // Omega < 1, dissolution kinetics
        {
            // Neutral mechanism: knu = 2.1E-09 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 2.1E-09 * exp(- 56.0 / R * (1.0/T - 1.0/298.15));

            // Acidic mechanism: k1 = 5.9E-06 * exp((-56000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.60)
            const auto kappa_acid = 5.9E-06 * exp(- 56.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.6);

            const auto kappa_diss = kappa_neu + kappa_acid;

            // Calculate the resulting mechanism function
            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa_diss * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;

        }
        else // Omega > 1, precipitation kinetics
        {
            // Interfacial energy in J m-2 # surf_energy = 0.06
            const auto surf_energy = 0.06;

            // Molecular volume in m3 # molvol = 3.4e-29
            const auto molecular_volume = 3.4e-29;

            // u = ( 16 * 3.14159 * surf_energy^3 * molvol ^ 2 ) / ( 3 * (1.38E-23 * TK)^3 ) 	# (after Fritz et al 2009)
            const auto u = 16 * 3.14159 * std::pow(surf_energy, 3.0) * std::pow(molecular_volume, 2.0) / (3 * std::pow(1.38E-23 * T.val, 3.0));

            // Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009) # J0 = 1E20
            const auto j0 = 1E20;

            // Critical saturation threshold # SRcrit = exp( sqrt( u / LOG( J0 ) ) )
            const auto Omega_crit = std::exp(sqrt(u / std::log10(j0)));

//            // Nucleation occurs only if saturation index above the threshold
//            if(Omega >= Omega_crit)
//            {
//                // Nucleation rate in nuclie sec-1 # nuclie = J0 * exp( -u / LOG(SRmin)^2 ) * time
//                const auto nuclie = j0 * exp(- u / std::pow(std::log10(Omega), 2.0));
//
//                if(nuclie * dt >= 1) // Condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
//                {
//                    // Number of growth units in critical nuclie radius # nj = ( 2 * u ) / (LOG(SRmin)^3) 		#
//                    const auto nj = 2 * u / std::pow(std::log10(Omega), 3.0);
//
//                    // Critical nuclie volume # vol = nj * molvol
//                    const auto vol = nj * molecular_volume;
//
//                    // Moles of nuclie formed # moles_nuc = ( -nuclie * vol ) / Mv
//                    res_nuc += - nuclie * vol / molar_volume;
//                }
//            }

            // Neutral mechanism: knu = 1.6E-11 * exp((-108000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 1.6E-11 * exp(- 108.0 / R * (1.0 / T - 1.0 / 298.15));

            // kpre = (-1) * knu
            const auto kappa_pre = - kappa_neu;

            // rate = S * m * Mm * kpre * (SRmin - 1)
            res_growth += f * ssa * nm * molar_mass * kappa_pre * std::abs(Omega - 1);

            // Precipitation kinetic includes kinetic growth and nucleation
            res += res_growth + res_nuc;

            // TODO: implement upper bound in precipitation kinetics
            // Do not precipitate more than the elements in solution ####
            // 320 maxMol = TOT("Fe(2)")
            // 330 IF (maxMol > TOT("C")) THEN maxMol = TOT("C")
            // 340 IF (maxMol < -moles) THEN moles = -maxMol

            // auto max_mol = tot_fe;
            // if(max_mol > tot_c) max_mol = tot_c;
            // if(max_mol < -res) res = - max_mol;
        }

        return res;
    };
    reaction_siderite.setRate(rate_func_siderite_olimse_dat);

    // Kaolinite: Al2Si2O5(OH)4
    // Al2Si2O5(OH)4 + 6.000H+ = 2.000Al+3 + 2.000H4SiO4 + 1.000H2O
    //     log_k     6.471
    //     delta_h -169.718    #kJ/mol        #01fia/nav
    //     -analytic -9.8589763E+2  -1.6937521E-1  5.5087963E+4  3.5699227E+2  -2.2447679E+6
    //     #References = LogK/DGf: Internal calculation; DHf/DHr: 01fia/nav; S°: 91rob/hem; Cp: 91rob/hem; V°: 95rob/hem;
    std::string eq_str_kaolinite = "Kaolinite + 6*H+ = 2*Al+++ + 2*SiO2(aq) + 5*H2O(l)";
    MineralReaction min_reaction_kaolinite = editor.addMineralReaction("Kaolinite")
            .setEquation(eq_str_kaolinite)
            .setSpecificSurfaceArea(11.8, "m2/g"); // S = 11.8 # average BET; suggested value in m2/g
    Reaction reaction_kaolinite = createReaction(min_reaction_kaolinite, system);
    reaction_kaolinite.setName("Kaolinite reaction");
    ReactionRateFunction rate_func_kaolinite_olimse_dat = [&min_reaction_kaolinite, &reaction_kaolinite, &system, &num_species, &R, &dt](const ChemicalProperties& properties) -> ChemicalScalar {

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0), res_growth(num_species, 0.0), res_nuc(num_species, 0.0);

        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_kaolinite.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_kaolinite.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        // const Index i_si = system.indexElementWithError("Si");
        // const auto tot_si = b(i_si);

        // Calculate TOT("Al")
        // const Index i_al = system.indexElementWithError("Al");
        // const auto tot_al = b(i_al);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_kaolinite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral].val;

        // Calculate activities for H+ and OH- species
        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        const Index i_oh = system.indexSpeciesWithError("OH-");
        double activity_oh = std::exp(lna(i_oh));

        // If (m <= 0) and (SRmin < 1) Then GoTo 350
        // If (SRmin = 1) Then GoTo 350
        if((nm <= 0 && Omega < 1) || Omega == 1) // the is no way to precipitate further
            return res;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area, average BET
        const auto ssa = min_reaction_kaolinite.specificSurfaceArea();

        // Molar volume # Mv = 9.876E-05, in m3/mol
        const auto molar_volume = 9.876E-05;

        if(Omega < 1) // Omega < 1, dissolution kinetics
        {
            // Neutral mechanism: knu = 1.1E-14 * exp((-38000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 1.1E-14 * exp(- 38.0 / R * (1.0/T - 1.0/298.15));

            // Acidic mechanism: k1 = 7.5E-12 * exp((-43000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H+") ^ 0.51)
            const auto kappa_acid = 7.5E-12 * exp(- 43.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.51);

            // OH- catalyzer mechanism (hydroxide promotion): k2 = 2.5E-11 * exp((-46000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.58)
            const auto kappa_oh = 2.5E-11 * exp(- 46.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_oh, 0.58);

            const auto kappa_diss = kappa_neu + kappa_acid + kappa_oh;

            // Calculate the resulting mechanism function # rate = S * m * Mm * k * ((1 - SRmin))
            res += f * ssa * nm * molar_mass * kappa_diss * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;

        }
        else // Omega > 1, precipitation kinetics
        {
            // Interfacial energy in J m-2 # surf_energy = 0.1
            const auto surf_energy = 0.1;

            // Thickness of one mineral layer in m # sheet_thick = 7e-10
            const auto sheet_thick = 7e-10;

            // Molecular volume in m3 # molvol = 1.64E-28
            const auto molecular_volume = 1.64E-28;

            // u = ( 2 * sqrt(3) * sheet_thick * surf_energy^2 * molvol ) / ( (1.38E-23 * TK)^2 ) 	# (after Fritz et al 2009)
            const auto u = 2 * std::sqrt(3) * sheet_thick * std::pow(surf_energy, 2.0) * molecular_volume / std::pow(1.38E-23 * T.val, 2.0);

            // Nucleation rate in nuclei kg-1 sec-1 (after Fritz et al 2009) # J0 = 1E20
            const auto j0 = 1E20;

            // Critical saturation threshold # SRcrit = exp( u / LOG( J0 ) )
            const auto Omega_crit = std::exp(u / std::log10(j0));

//            // Nucleation occurs only if saturation index above the threshold
//            if(Omega >= Omega_crit)
//            {
//                // Nucleation rate in nuclie sec-1 # nuclie = J0 * exp( -u / LOG(SRmin) ) * Time
//                const auto nuclie = j0 * exp(- u / std::log10(Omega));
//
//                if(nuclie * dt >= 1) // Condition that rate needs to be bigger than 1 nucl per sec (Arbitrary)
//                {
//                    // Number of growth units in critical nuclie radius # 250 nj = u / (LOG(SRmin)^2)
//                    const auto nj = 2 * u / std::pow(std::log10(Omega), 3.0);
//
//                    // Critical nuclie volume # vol = nj * molvol
//                    const auto vol = nj * molecular_volume;
//
//                    // Moles of nuclie formed # moles_nuc = ( -nuclie * vol ) / Mv
//                    res_nuc += - nuclie * vol / molar_volume;
//                }
//            }

            // knu = 5.5E-13 * exp((-66000 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 5.5E-13 * exp(- 66.0 / R * (1.0 / T - 1.0 / 298.15));

            // kpre = (-1) * knu
            const auto kappa_pre = - kappa_neu;

            // moles_growth = S * m * Mm * kpre * (ABS((SRmin ^ 0.06) - 1) ^ 1.68) * Time
            res_growth += f * ssa * nm * molar_mass * kappa_pre * std::pow(std::abs(std::pow(Omega, 0.06) - 1), 1.68);

            res += res_growth + res_nuc;

            // TODO: implement upper bound in precipitation kinetics
            // Do not precipitate more than the elements in solution ####
            // 360 maxMol = TOT("Al")/2
            // 370 IF (maxMol > TOT("Si")/2) THEN maxMol = TOT("Si")/2
            // 380 IF (maxMol < -moles) THEN moles = -maxMol

            // auto max_mol = tot_al / 2;
            // if(max_mol > tot_si / 2) max_mol = tot_si / 2;
            // if(max_mol < -res) res = - max_mol;
        }

        return res;
    };
    reaction_kaolinite.setRate(rate_func_kaolinite_olimse_dat);

    // Quartz(alpha)
    //SiO2 + 2.000H2O = 1.000H4SiO4
    //     log_k    -3.737
    //     delta_h  21.166     #kJ/mol        #82ric/bot
    //     -analytic -7.5895338E+1  -1.5422139E-2  1.5615589E+3  2.9087273E+1  -4.0514987E+4
    //     #References = LogK/DGf: Internal calculation; DHf/DHr: 82ric/bot; S°: 82ric/bot; Cp: 82ric/bot; V°: 95rob/hem;
    std::string eq_str_quartz = "Quartz + H2O(l) = H+ + HSiO3-";
    //std::string eq_str_quartz = "Quartz = SiO2(aq)";
    MineralReaction min_reaction_quartz = editor.addMineralReaction("Quartz")
            .setEquation(eq_str_quartz)
            .setSpecificSurfaceArea(0.03, "m2/g"); // S = 0.03 # average BET; suggested value in m2/g
    Reaction reaction_quartz = createReaction(min_reaction_quartz, system);
    reaction_quartz.setName("Quartz reaction");
    ReactionRateFunction rate_func_quartz_olimse_dat = [&min_reaction_quartz, &reaction_quartz, &system, &num_species, &R, &dt](const ChemicalProperties& properties) -> ChemicalScalar {

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);

        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Calculate the saturation index of the mineral
        const auto lnK = reaction_quartz.lnEquilibriumConstant(properties);
        const auto lnQ = reaction_quartz.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        auto Omega = exp(lnOmega).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0 = reaction_quartz.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        // const Index i_si = system.indexElementWithError("Si");
        // const auto tot_si = b(i_si);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_quartz.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral].val;
        auto nm0 = n0[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_oh = system.indexSpeciesWithError("OH-");
        double activity_oh = std::exp(lna(i_oh));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // If (m <= 0) and (SRmin < 1) Then GoTo 240
        // If (SRmin = 1) Then GoTo 240
        if((nm <= 0 && Omega < 1) || (Omega == 1)) // the is no way to precipitate further
            return res;

        // Specific surface area, average BET
        const auto ssa = min_reaction_quartz.specificSurfaceArea();

        if(Omega < 1) // Omega < 1, dissolution kinetics
        {
            // Neutral mechanism: knu = 6.42E-14 * exp((-76700 / 8.314) * ((1 / TK) - (1 / 298.15)))
            const auto kappa_neu = 6.42E-14 * exp(- 76.7 / R * (1.0/T - 1.0/298.15));

            // OH- catalyzer mechanism (hydroxide promotion): k1 = 0.000000000192 * exp((-80000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("OH-") ^ 0.339)
            const auto kappa_oh = 0.000000000192 * exp(- 80.0 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_oh, 0.339);

            const auto kappa_diss = kappa_neu + kappa_oh;

            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - SRmin) # by default
            res += f * ssa * nm * molar_mass * pow(nm / nm0, 2.0 / 3.0) * kappa_diss * (1 - Omega);

            // TODO: Do not dissolve more than what is available
            // IF (moles > M) THEN moles = M
            // if (res > nm.val) res = nm.val;
        }
        else // Omega > 1, precipitation kinetics
        {
            // If (m <= 1e-5) then GoTo 170
            if(nm > 1e-5) // start precipitation
            {
                // knu = 3.24e-12
                // kpre = (-1) * knu
                const auto kappa_pre = - 3.24e-12;

                // rate = S * m * Mm * kpre * (ABS((SRmin ^ 4.58) - 1) ^ 0.54)
                res += f * ssa * nm * molar_mass * kappa_pre * std::pow(std::abs(std::pow(Omega, 4.58) - 1), 0.54);
            }
            else // start nucleation
            {
                // Set nucleation rate
                res += - 1e-10 * f;
            }

            // TODO: implement upper bound in precipitation kinetics
            // auto max_mol = tot_si;
            // if(max_mol < -res) res = - max_mol;
        }
        return res;

    };
    reaction_quartz.setRate(rate_func_quartz_olimse_dat);

    // Create the ReactionSystem instances
    ReactionSystem reactions(system, {reaction_calcite, reaction_siderite, reaction_daphnite, reaction_kaolinite, reaction_quartz});
    //ReactionSystem reactions(system, {reaction_calcite, reaction_siderite, reaction_daphnite, reaction_kaolinite, reaction_quartz});

    Partition partition(system);
    partition.setKineticSpecies(std::vector<std::string>{"Calcite", "Siderite", "Daphnite,14A", "Kaolinite", "Quartz"});

    double water_kg = 1.00;
    EquilibriumInverseProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    // Seawater composition
    problem_ic.add("H2O", water_kg, "kg");
    problem_ic.add("HCO3-", 1441.11, "mg"); // 1441.11 mg/l
    problem_ic.add("Ca++", 310, "mg"); // 310 mg/l
    problem_ic.add("Cl-", 19100, "mg"); // 19100 charge
    problem_ic.add("Fe++", 0.0034, "mg"); // 0.0034 mg/l
    problem_ic.add("Mg++", 1240, "mg"); // 1240 mg/l
    problem_ic.add("Na+", 11400, "mg"); // 11400 mg/l
    problem_ic.add("SO4--", 2420, "mg"); // 2420 mg/l
    // Sulfide spike
    problem_ic.add("S", 100, "mg"); // 100 mg/l
    // Equilibrium phases
    problem_ic.add("Calcite", 0.0, "mol"); // Calcite 0 0
    problem_ic.add("Pyrrhotite", 0.0, "mol"); // Pyrrhotite 0 0
    problem_ic.pH(7.1015);
    problem_ic.pE(-3.6887);

    ChemicalState state_ic = equilibrate(problem_ic);

    state_ic.setSpeciesAmount("Calcite", 0.1, "mol"); // MM(100.09) = 713.5 g / mol
    state_ic.setSpeciesAmount("Daphnite,14A", 0.1, "mol"); // MM(Daphnite) = 713.5 g / mol
    state_ic.setSpeciesAmount("Siderite", 0.1, "mol"); // MM(Siderite) = 115.86 g / mol
    state_ic.setSpeciesAmount("Kaolinite", 0.1, "mol"); // MM(Kaolinite) = 258.071 g / mol
    state_ic.setSpeciesAmount("Quartz", 1.0, "mol"); // MM(Quartz) = 60.083 g / mol

    // Set initial value of minerals reacting
    reaction_calcite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_daphnite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_siderite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_kaolinite.setInitialAmounts(state_ic.speciesAmounts());
    reaction_quartz.setInitialAmounts(state_ic.speciesAmounts());

    // Initialize kinetic path
    KineticPath path(reactions, partition);
    path.setOptions(kineticpath_options);

    // -------------------------------------------------------------------------------------------------//
    // Approach III: solve the path with embedded in the KineticPath uniform time-stepping procedure
    // -------------------------------------------------------------------------------------------------//
    ChemicalOutput output = path.output();
    output.add("time(units=s)");
    output.add("pH");
    output.add("speciesAmount(Ca++)");
    output.add("speciesAmount(Mg++)");
    output.add("speciesAmount(Na+)");
    output.add("speciesAmount(Fe++)");
    output.add("speciesAmount(Fe+++)");
    output.add("speciesAmount(SO4--)");
    output.add("speciesAmount(HS-)");
    output.add("speciesAmount(HCO3-)");
    output.add("speciesAmount(H+");
    output.add("speciesAmount(Al+++");
    output.add("speciesAmount(SiO2(aq))");
    output.add("speciesAmount(HSiO3-");
    output.add("speciesAmount(Pyrrhotite)");
    output.add("speciesAmount(Calcite)");
    output.add("speciesAmount(Quartz)");
    output.add("speciesAmount(Kaolinite)");
    output.add("speciesAmount(Daphnite,14A)");
    output.add("speciesAmount(Siderite)");
    output.filename(filename);

    tic(TRANSPORT)

    path.solve(state_ic, params.t0, params.dt, params.n, "second");

    double total_time = toc(TRANSPORT);
    std::cout << filename << std::endl;
    std::cout << "total_time = " << total_time << std::endl;

}

