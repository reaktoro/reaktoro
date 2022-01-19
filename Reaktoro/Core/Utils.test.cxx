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
#include <Reaktoro/Core/Utils.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing CoreUtils", "[CoreUtils]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto mmH2O = system.species().get("H2O(aq)").molarMass(); // kg/mol
    const auto iH2O = system.species().index("H2O(aq)");

    real actual = {};
    real expected = {};

    SECTION("Testing computeSpeciesAmount")
    {
        actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "mol");
        expected = 1.0;
        CHECK( actual == expected );

        actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "mmol");
        expected = 0.001;
        CHECK( actual == expected );

        actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "kg");
        expected = 1.0 / mmH2O;
        CHECK( actual == Approx(expected) );

        actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "g");
        expected = 0.001 / mmH2O;
        CHECK( actual == Approx(expected) );
    }

    auto species = system.species();
    auto elements = system.elements();

    auto components_num = elements.size() + 1;

    SECTION("Testing assembleFormulaMatrix for the full system")
    {
        auto A = detail::assembleFormulaMatrix(species, elements);

        CHECK( A.cols() == species.size() );
        CHECK( A.rows() == components_num );

        MatrixXd A_system = system.formulaMatrix();

        for(Index i = 0; i < A.rows(); i++)
            CHECK(A.row(i) == A_system.row(i));
    }

    SECTION("Testing assembleFormulaMatrix for the aqueous phase")
    {
        SpeciesList species_aqueous = species.withAggregateState(AggregateState::Aqueous);

        auto A_aq = detail::assembleFormulaMatrix(species_aqueous, elements);

        CHECK( A_aq.cols() == species_aqueous.size() );
        CHECK( A_aq.rows() == components_num         );

        MatrixXd A_aq_expected {
            {2,  1,  1,  2,  0,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0},
            {1,  0,  1,  0,  2,  0,  0,  0,  0,  1,  0,  0,  2,  3,  3,  0,  0,  2,  0},
            {0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0},
            {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0},
            {0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0},
            {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0},
            {0,  1, -1,  0,  0,  1, -1,  0,  0,  0,  2,  2,  0, -1, -2,  0,  0,  0, -1}
        };

        for(Index i = 0; i < A_aq.rows(); i++)
            CHECK(A_aq.row(i) == A_aq_expected.row(i));
    }

    SECTION("Testing assembleFormulaMatrix for the gaseous phase")
    {
        SpeciesList species_gaseous = species.withAggregateState(AggregateState::Gas);

        auto A_gas = detail::assembleFormulaMatrix(species_gaseous, elements);

        CHECK( A_gas.cols() == species_gaseous.size() );
        CHECK( A_gas.rows() == components_num         );

        MatrixXd A_gas_expected {
            {0,  0,  2,  2,  4,  0},
            {1,  0,  0,  0,  1,  1},
            {2,  2,  0,  1,  0,  1},
            {0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0},
            {0,  0,  0,  0,  0,  0},
        };

        for(Index i = 0; i < A_gas.rows(); i++)
            CHECK(A_gas.row(i) == A_gas_expected.row(i));
    }

    SECTION("Testing assembleFormulaMatrix for the gaseous phase")
    {
        SpeciesList species_minerals = species.withAggregateState(AggregateState::Solid);

        auto A_min = detail::assembleFormulaMatrix(species_minerals, elements);

        CHECK( A_min.cols() == species_minerals.size() );
        CHECK( A_min.rows() == components_num          );

        MatrixXd A_min_expected {
            {0,  0,  0,  0},
            {0,  1,  1,  0},
            {0,  3,  3,  2},
            {1,  0,  0,  0},
            {0,  0,  1,  0},
            {0,  0,  0,  1},
            {1,  0,  0,  0},
            {0,  1,  0,  0},
            {0,  0,  0,  0},
        };

        for(Index i = 0; i < A_min.rows(); i++)
            CHECK(A_min.row(i) == A_min_expected.row(i));
    }

}
