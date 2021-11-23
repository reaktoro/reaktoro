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
//   • Svetlana Kyas (23 November 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
#include <Reaktoro/Core/Utils.hpp>

using namespace Reaktoro;

int main()
{
    // Initialize the database
    auto dbphreeqc = PhreeqcDatabase("phreeqc.dat");

    for(auto s : dbphreeqc.species()) {
        std::cout << s.name() << std::endl;
    }
    // Define ion exchange species list
    // Expected species: Hfo_psi Hfo_psib Hfo_psid Hfo_sOH Hfo_wOH Hfo_wF
    // Hfo_wH2BO3 Hfo_wH2PO4 Hfo_wH3SiO4 Hfo_wHCO3 Hfo_wOFeOH Hfo_wOSrOH
    SpeciesList slist = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed).withCharge(0.0);

    for(auto s : slist) {
        std::cout << s.name() << " ";
    }
    std::cout << std::endl;

    // Define an surface complexation phase
    SurfaceComplexationPhase complexation_phase(detail::extractNames(slist));
    complexation_phase.setActivityModel(ActivityModelSurfaceComplexationDDL());

    SurfaceComplexationPhase complexation_phase_Hfo("Hfo_sOH Hfo_wOH");
    complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationDDL());

    for(auto s : complexation_phase_Hfo.species()) {
        std::cout << s << " ";
    }
    std::cout << std::endl;

    return 0;
}