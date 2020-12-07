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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);
    StringList selected_elements = "H Cl S O Ba Ca Sr Na K Mg C Si";

    editor.addAqueousPhaseWithElements(selected_elements)
            .setChemicalModelDebyeHuckel(dhModel);
    editor.addMineralPhase("Barite");

    // Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Create reaction for barite mineral
    std::string eq_str_barite = "Barite = SO4-- + Ba++";
    MineralReaction min_reaction_barite = editor.addMineralReaction("Barite")
            .setEquation(eq_str_barite)
            .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol")
            .setSpecificSurfaceArea(0.006, "m2/g");

    // Create reaction for barite
    Reaction reaction_barite = createReaction(min_reaction_barite, system);
    reaction_barite.setName("Barite kinetics");
    // Ionic strength function
    const auto I = ChemicalProperty::ionicStrength(system);
    // Barite kinetic rate function
    ReactionRateFunction rate_func_barite_shell = [&min_reaction_barite, &reaction_barite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate using specified surface area
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
        const auto Omega = exp(lnOmega);

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

        // The number of moles of the mineral
        auto nm = n[imineral];
        auto nm0 = n0[imineral];

        // Prevent negative mole numbers here for the solution of the ODEs
        nm.val = std::max(nm.val, 0.0);
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
            res = 0;
        // S = 0.006 # average BET from 16zhe/did ; suggested value in m2/g
        const auto ssa = 0.006;

        // If (SRmin > 1) Then GoTo 130
        if(Omega > 1) // precipitation kinetics
        {
            /*
             * ########## start precipitation bloc ##########
            130 If (m <= 1e-5) then GoTo 170
            140 knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
            150 kpre = (-1) * knu
            160 rate = S * m * Mm * kpre * (SRmin - 1)
            170 GoTo 190
            #start nucleation
            180 rate = -1e-12
            190 moles = rate * Time
            200 REM Do not precipitate more than the elements in solution
            210 maxMol = TOT("Ba")
            220 IF (maxMol > TOT("S(6)")) THEN maxMol = TOT("S(6)")
            230 IF (maxMol < -moles) THEN moles = -maxMol
            ########## end precipitation bloc ##########
             */

            // If (m <= 1e-5) then GoTo 170
            if(nm > 1e-5)
            {
                // Calculate the rate constant
                // knu = 2.18E-09 * exp((-22000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (10 ^ (0.6 * MU^0.5))
                // kpre = (-1) * knu
                const auto kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15))
                                       * pow(10, 0.6 * pow(ionic_strength, 0.5));
                //const auto kappa_pre = - 2.18E-09 * exp(- 22 / R * (1.0/T - 1.0/298.15));

                // rate = S * m * Mm * kpre * (SRmin - 1)
                res += f * ssa * nm.val * molar_mass * kappa_pre * (Omega - 1);
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
            /*
             * ########## start dissolution bloc ##########
            70 k1 = 2.75E-08 * exp((-25000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.03) * (10 ^ (0.6 * MU^0.5))
            80 k = k1
            90 rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - (SRmin ^ 0.2))
            100 moles = rate * Time
            110 IF (moles > M) THEN moles = M
            120 GoTo 240
            ########## end dissolution bloc ##########
             */

            // k1 = 2.75E-08 * exp((-25000 / 8.314) * ((1 / TK) - (1 / 298.15))) * (ACT("H3O+") ^ 0.03) * (10 ^ (0.6 * MU^0.5))
            const auto kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15))
                               * std::pow(activity_h, 0.03)
                               * pow(10, 0.6 * pow(ionic_strength, 0.5));

            // rate = S * m * Mm * ((m/m0)^(2/3)) * k * (1 - (SRmin ^ 0.2))
            res += f * ssa * nm.val * molar_mass * pow(nm.val / nm0, 2 / 3) * kappa * (1 - pow(Omega, 0.2));

        }

        return res;

    };
    reaction_barite.setRate(rate_func_barite_shell);

    // Create the ReactionSystem instances
    ReactionSystem reactions(system, {reaction_barite});
    //ReactionSystem reactions(editor); // will fetch predifined reaction rates

    Partition partition(system);
    partition.setKineticPhases({"Barite"});

    double water_kg = 1.00;
    double T = 60;
    double P = 200;

    EquilibriumInverseProblem problem_ic_fw(system);
    problem_ic_fw.setTemperature(T, "celsius");
    problem_ic_fw.setPressure(P, "atm");
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

    ChemicalState state_ic = equilibrate(problem_ic_fw);
    state_ic.setSpeciesAmount("Barite", 0.1, "mcmol");
    state_ic.scaleVolume(1.0, "m3");

    // Set initial value of Barite
    reaction_barite.setInitialAmounts(state_ic.speciesAmounts());

    // Define the first boundary condition (with seawater)
    EquilibriumInverseProblem problem_bc_sw(system);
    problem_bc_sw.setTemperature(T, "celsius");
    problem_bc_sw.setPressure(P, "atm");
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
    ChemicalState state_bc = equilibrate(problem_bc_sw);
    state_bc.scaleVolume(1.0, "m3");

    std::cout << state_bc << std::endl;

    KineticPath path(reactions, partition);

//    // -------------------------------------------------------------------------------------------------//
//    // Approach I: solve the path on one interval adn let CVODE treat interval [t0, t0 + dt] adaptively
//    // -------------------------------------------------------------------------------------------------//
//
//    ChemicalOutput output = path.output();
//    output.add("time(units=s)");
//    output.add("speciesAmount(Barite units=mol)", "Barite");
//    output.add("speciesAmount(Ba++ units=mol)", "Ba++");
//    output.add("speciesAmount(SO4-- units=mol)", "SO4--");
//    output.filename("kineticpath-barite-precipitation-one-interval.txt");
//
//    double t0 = 0;
//    double dt = 50.0;
//    state_ic = state_ic + state_bc;
//    path.solve(state_ic, t0, t0 + dt, "days");


//    // -------------------------------------------------------------------------------------------------//
//    // Approach II: solve the path on one interval adn let CVODE treat interval [t0, t0 + dt] adaptively
//    // -------------------------------------------------------------------------------------------------//
//
//    state_ic = state_ic + state_bc;
//
//    double t0 = 0;
//    double dt = 5.0;
//    double tfinal = 50.0;
//    int n = tfinal / dt;
//
//    for(int i = 0; i < n; i++){
//
//        ChemicalOutput output = path.output();
//        output.add("time(units=s)");
//        output.add("speciesAmount(Barite units=mol)", "Barite");
//        output.add("speciesAmount(Ba++ units=mol)", "Ba++");
//        output.add("speciesAmount(SO4-- units=mol)", "SO4--");
//        output.filename("kineticpath-barite-precipitation-" + std::to_string(i) + ".txt");
//
//        path.solve(state_ic, t0, t0 + dt, "days");
//        t0 += dt;
//    }

    // -------------------------------------------------------------------------------------------------//
    // Approach III: solve the path with embedded in the KineticPath uniform time-stepping procedure
    // -------------------------------------------------------------------------------------------------//
    ChemicalOutput output = path.output();
    output.add("time(units=s)");
    output.add("speciesAmount(Barite units=mol)", "Barite");
    output.add("speciesAmount(Ba++ units=mol)", "Ba++");
    output.add("speciesAmount(SO4-- units=mol)", "SO4--");
    output.filename("kineticpath-barite-kinetics-n-steps.txt");

    double t0 = 0;
    double dt = 5.0;
    double tfinal = 50.0;
    int n = tfinal / dt;
    state_ic = state_ic + state_bc;
    path.solve(state_ic, t0, dt, n, "day");

}
