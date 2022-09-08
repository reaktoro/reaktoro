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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumDims", "[EquilibriumDims]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto Ne = system.elements().size();
    const auto Nn = system.species().size();

    EquilibriumSpecs specs(system);

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );  // the number of elements in the chemical system.
        CHECK( dims.Nn == Nn     );  // the number of species in the chemical system.
        CHECK( dims.Np == 0      );  // the number of *p* control variables (temperature, pressure, amounts of explicit titrants, and custom variables).
        CHECK( dims.Nq == 0      );  // the number of *q* control variables (amounts of implicit titrants).
        CHECK( dims.Nv == 0      );  // the number of equations constraints in the chemical equilibrium problem.
        CHECK( dims.Nr == 0      );  // the number of reactivity constraints (i.e., *restricted reactions*) in the chemical equilibrium problem.
        CHECK( dims.Nb == Ne + 1 ); // elements, charge  // the number of components (electric charge, chemical elements, extent of restricted reactions) in the chemical equilibrium problem (equivalent to `Ne + 1 + Nr`).
        CHECK( dims.Nt == 0      );  // the number of substances for which the chemical system is open to (the number of explicit and implicit titrants).
        CHECK( dims.Nx == Nn     );  // the number of variables *x* in *x = (n, q)* (equivalent to `Nn + Nq`).
        CHECK( dims.Nu == Nn     );  // the number of unknown variables in the chemical equilibrium problem (equivalent to `Nn + Np + Nq`).
        CHECK( dims.Nw == 2      );  // the number of input variables *w* in the chemical equilibrium problem.
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 1      ); // P
        CHECK( dims.Nq == 0      );
        CHECK( dims.Nv == 1      ); // V
        CHECK( dims.Nr == 0      );
        CHECK( dims.Nb == Ne + 1 ); // elements, charge
        CHECK( dims.Nt == 0      );
        CHECK( dims.Nx == Nn     );
        CHECK( dims.Nu == Nn + 1 );
        CHECK( dims.Nw == 2      ); // T, V
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 2      ); // T, P
        CHECK( dims.Nq == 0      );
        CHECK( dims.Nv == 2      ); // V, U
        CHECK( dims.Nr == 0      );
        CHECK( dims.Nb == Ne + 1 ); // elements, charge
        CHECK( dims.Nt == 0      );
        CHECK( dims.Nx == Nn     );
        CHECK( dims.Nu == Nn + 2 );
        CHECK( dims.Nw == 2      ); // V, U
    }

    WHEN("temperature, pressure, and pH are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 0      );
        CHECK( dims.Nq == 1      ); // [H+]
        CHECK( dims.Nv == 0      );
        CHECK( dims.Nr == 0      );
        CHECK( dims.Nb == Ne + 1 ); // elements, charge
        CHECK( dims.Nt == 1      ); // [H+]
        CHECK( dims.Nx == Nn + 1 );
        CHECK( dims.Nu == Nn + 1 );
        CHECK( dims.Nw == 3      ); // T, P, pH
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input variables")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 2      );
        CHECK( dims.Nq == 1      ); // [CO2]
        CHECK( dims.Nv == 2      ); // V, S
        CHECK( dims.Nr == 0      );
        CHECK( dims.Nb == Ne + 1 ); // elements, charge
        CHECK( dims.Nt == 1      ); // [CO2]
        CHECK( dims.Nx == Nn + 1 );
        CHECK( dims.Nu == Nn + 3 );
        CHECK( dims.Nw == 3      ); // T, P, a(CO2)
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 2      ); // [CO2], [CH4]
        CHECK( dims.Nq == 2      ); // [H+], [e-]
        CHECK( dims.Nv == 2      ); // V, U
        CHECK( dims.Nr == 0      );
        CHECK( dims.Nb == Ne + 1 ); // elements, charge
        CHECK( dims.Nt == 4      );
        CHECK( dims.Nx == Nn + 2 );
        CHECK( dims.Nu == Nn + 4 );
        CHECK( dims.Nw == 6      ); // T, P, V, U, pH, pE
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables, there are reactivity constraints")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");
        specs.addReactivityConstraint({ "xi1", VectorXd::Random(Nn), {} });
        specs.addReactivityConstraint({ "xi2", VectorXd::Random(Nn), {} });

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne     );
        CHECK( dims.Nn == Nn     );
        CHECK( dims.Np == 2      ); // [CO2], [CH4]
        CHECK( dims.Nq == 2      ); // [H+], [e-]
        CHECK( dims.Nv == 2      ); // V, U
        CHECK( dims.Nr == 2      );
        CHECK( dims.Nb == Ne + 3 ); // elements, charge, xi1, xi2
        CHECK( dims.Nt == 4      );
        CHECK( dims.Nx == Nn + 2 );
        CHECK( dims.Nu == Nn + 4 );
        CHECK( dims.Nw == 6      ); // T, P, V, U, pH, pE
    }
}
