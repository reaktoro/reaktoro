// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Svetlana Kyas (30 August 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

const auto T = 25.0 + 273.15; // temperature in K
const auto P = 1.0 * 1e5;    // pressure in Pa

auto speciesListToStringList(const SpeciesList& specieslist) -> StringList
{
    std::vector<std::string> speciesvector;
    for (auto species : specieslist)
        speciesvector.push_back(species.name());

    return StringList(speciesvector);
}
int main()
{
    // Initialize Phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define aqueous phase
    AqueousPhase aqueousphase(speciate("H O C Ca Na Mg Cl"));
    aqueousphase.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond("CO2")
            ));

    // Fetch species for ion-exchange modeling
    SpeciesList exchange_species = db.species().withAggregateState(AggregateState::IonExchange);
    // Separate exchangres (e.g., X-) from exchange species species (NaX, KX, NaY, NaY, ect)
    SpeciesList exchanger_species = filter(exchange_species, [](const Species& s){ return std::abs(s.charge());});
    SpeciesList ionexchange_species = filter(exchange_species, [](const Species& s){ return !s.charge();});

    // The exchanger (exchanging site) phase X with exchange species X-
    GenericPhase exchanger_phase(speciesListToStringList(exchanger_species));
    exchanger_phase.setName("ExchangerPhase");
    exchanger_phase.setStateOfMatter(StateOfMatter::Solid);
    exchanger_phase.setAggregateState(AggregateState::IonExchange);
    exchanger_phase.setActivityModel(ActivityModelIdealSolution());
    exchanger_phase.setIdealActivityModel(ActivityModelIdealSolution());

    // Ion exchange phase (the cation exchange complex), containing species NaX, KX, CaX2 etc
    GenericPhase ionexchange_phase(speciesListToStringList(ionexchange_species));
    ionexchange_phase.setName("IonExchangePhase");
    ionexchange_phase.setStateOfMatter(StateOfMatter::Solid);
    ionexchange_phase.setAggregateState(AggregateState::IonExchange);
    ionexchange_phase.setActivityModel(ActivityModelIonExchangeGainesThomas());
    ionexchange_phase.setIdealActivityModel(ActivityModelIdealSolution());

    // Construct the chemical system
    ChemicalSystem system(db, aqueousphase, exchanger_phase, ionexchange_phase);

    // Specify conditions to be satisfied at the chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.charge();
    specs.openTo("Cl-");

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(60.0, "celsius");
    conditions.pressure(100.0, "bar");
    conditions.charge(0, "mol");

    // Define initial equilibrium state
    ChemicalState solutionstate(system);
    solutionstate.setTemperature(T, "celsius");
    solutionstate.setPressure(P, "bar");
    solutionstate.setSpeciesMass("H2O"    , 1.00, "kg");
    solutionstate.setSpeciesAmount("Na+"  , 1.10, "mol");
    solutionstate.setSpeciesAmount("Mg+2" , 0.48, "mol");
    solutionstate.setSpeciesAmount("Ca+2" , 1.90, "mol");
    solutionstate.setSpeciesAmount("X-"   , 0.06, "mol");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);
    solver.solve(solutionstate, conditions);

    // Output the chemical state to a text file
    solutionstate.output("state.txt");

    return 0;
}
