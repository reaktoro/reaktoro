// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystemAux.hpp>
using namespace Reaktoro;

int main()
{
    // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT

    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)");
    solution.setActivityModel(ActivityModelPitzerHMW());

    MineralPhase calcite("Calcite");

    Reaction reaction1 = db.reaction("H2O(aq) = H+ + OH-");
    Reaction reaction2 = db.reaction("Calcite = Ca+2 + CO3-2");

    Params params = Params::embedded("PalandriKharaka.yaml");

    MineralReactions mineralreactions("Calcite");
    mineralreactions.setRateModel(ReactionRateModelPalandriKharaka(params));

    // ChemicalSystem system(db, solution, calcite);
    ChemicalSystem system(db, mineralreactions, solution, calcite, reaction1, reaction2);

    ChemicalState state(system);
    state.temperature(50.0, "celsius");
    state.pressure(2.0, "bar");
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Calcite", 1.0, "mol");

    std::cout << state << std::endl;

    for(auto reaction : system.reactions())
        std::cout << reaction.equation() << std::endl;

    for(auto phase : system.phases())
        std::cout << phase.name() << std::endl;
}
