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

    const auto Ne = system.elements().size() + 1;
    const auto Nn = system.species().size();

    EquilibriumSpecs specs(system);

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne ); // number of elements in the chemical system
        CHECK( dims.Nb == Ne ); // number of components in the chemical equilibrium problem
        CHECK( dims.Nn == Nn ); // number of species in the chemical system
        CHECK( dims.Np == 0 );  // number of *p* control variables (temperature, pressure, and amounts of explicit titrants when these are introduced unknowns)
        CHECK( dims.Nq == 0 );  // number of *q* control variables (the amounts of the implicit titrants when these are introduced unknowns)
        CHECK( dims.Nt == 0 );  // number of substances for which the chemical system is open to
        CHECK( dims.Nx == Nn ); // number of variables *x* in *x = (n, q)* (equivalent to `Nn + Nq`)
        CHECK( dims.Nu == Nn ); // number of unknown variables in the chemical equilibrium problem (equivalent to `Nn + Np + Nq`)
        CHECK( dims.Nw == 2  ); // number of input variables in the chemical equilibrium problem.
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne );
        CHECK( dims.Nb == Ne );
        CHECK( dims.Nn == Nn );
        CHECK( dims.Np == 1  ); // P
        CHECK( dims.Nq == 0  );
        CHECK( dims.Nt == 0  );
        CHECK( dims.Nx == Nn );
        CHECK( dims.Nu == Nn + 1 );
        CHECK( dims.Nw == 2  ); // T, V
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne );
        CHECK( dims.Nb == Ne );
        CHECK( dims.Nn == Nn );
        CHECK( dims.Np == 2  ); // T, P
        CHECK( dims.Nq == 0  );
        CHECK( dims.Nt == 0  );
        CHECK( dims.Nx == Nn );
        CHECK( dims.Nu == Nn + 2 );
        CHECK( dims.Nw == 2  ); // V, U
    }

    WHEN("temperature, pressure, and pH are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne );
        CHECK( dims.Nb == Ne );
        CHECK( dims.Nn == Nn );
        CHECK( dims.Np == 0 );
        CHECK( dims.Nq == 1 ); // [H+]
        CHECK( dims.Nt == 1 ); // [H+]
        CHECK( dims.Nx == Nn + 1 );
        CHECK( dims.Nu == Nn + 1 );
        CHECK( dims.Nw == 3  ); // T, P, pH
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input variables")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        EquilibriumDims dims(specs);

        CHECK( dims.Ne == Ne );
        CHECK( dims.Nb == Ne );
        CHECK( dims.Nn == Nn );
        CHECK( dims.Np == 2  );
        CHECK( dims.Nq == 1  ); // [CO2]
        CHECK( dims.Nt == 1  ); // [CO2]
        CHECK( dims.Nx == Nn + 1 );
        CHECK( dims.Nu == Nn + 3 );
        CHECK( dims.Nw == 3  ); // T, P, a(CO2)
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

        CHECK( dims.Ne == Ne );
        CHECK( dims.Nb == Ne );
        CHECK( dims.Nn == Nn );
        CHECK( dims.Np == 2  ); // [CO2], [CH4]
        CHECK( dims.Nq == 2  ); // [H+], [e-]
        CHECK( dims.Nt == 4  );
        CHECK( dims.Nx == Nn + 2 );
        CHECK( dims.Nu == Nn + 4 );
        CHECK( dims.Nw == 6 ); // T, P, V, U, pH, pE
    }
}
