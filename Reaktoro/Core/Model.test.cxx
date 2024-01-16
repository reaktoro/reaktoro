// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/Model.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Model class", "[Model]")
{
    auto K = 3.0;

    Data params;
    params["K"] = K;

    SECTION("Using ModelEvaluator")
    {
        auto evalfn = [=](real& res, real x, real y)
        {
            res = K*x*y;
        };

        Model<real(real, real)> model(evalfn, params);

        CHECK( model.initialized() );
        CHECK( model.evaluatorFn() );
        CHECK( model.calculatorFn() );

        CHECK( model.params().isDict() );
        CHECK( model.params()["K"].asFloat() == K );

        const auto x = 3.0;
        const auto y = 7.0;

        CHECK( model(x, y) == Approx(3.0 * x * y) );
    }

    SECTION("Using ModelCalculator")
    {
        auto calcfn = [=](real x, real y)
        {
            return K*x*y;
        };

        Model<real(real, real)> model(calcfn, params);

        CHECK( model.initialized() );
        CHECK( model.evaluatorFn() );
        CHECK( model.calculatorFn() );

        CHECK( model.params().isDict() );
        CHECK( model.params()["K"].asFloat() == K );

        const auto x = 3.0;
        const auto y = 7.0;

        CHECK( model(x, y) == Approx(3.0 * x * y) );
    }

    SECTION("Using ModelCalculator with params")
    {
        auto A = 1.0;
        auto B = 2.0;

        Data params;
        params["A"] = A;
        params["B"] = B;

        auto calcfn = [=](real x, real y) -> real
        {
            return A*x + B*y;
        };

        auto model = Model<real(real, real)>(calcfn, params);

        CHECK( model.params().isDict() );
        CHECK( model.params()["A"].asFloat() == A );
        CHECK( model.params()["B"].asFloat() == B );
    }

    SECTION("Using Model::Constant")
    {
        auto model = Model<real(real, real)>::Constant("MyConstantModel", K);

        CHECK( model.initialized() );
        CHECK( model.evaluatorFn() );
        CHECK( model.calculatorFn() );

        CHECK( model.params().isDict() );
        CHECK( model.params()["MyConstantModel"]["Value"].asFloat() == K );

        const auto x = 3.0;
        const auto y = 7.0;

        CHECK( model(x, y) == Approx(3.0) );
    }
}
