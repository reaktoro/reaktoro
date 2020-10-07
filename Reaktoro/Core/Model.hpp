// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

template<typename Signature>
class Model;

/// The functional signature of functions that evaluates properties.
template<typename ResultRef, typename... Args>
using ModelEvaluator = Fn<void(ResultRef res, Args... args)>;

/// The functional signature of functions that calculates properties.
template<typename Result, typename... Args>
using ModelCalculator = Fn<Result(Args... args)>;

/// The functional signature of functions that creates ModelEvaluator function objects.
template<typename ResultRef, typename... Args>
using ModelCreator = Fn<ModelEvaluator<ResultRef, Args...>(const Params& params)>;

/// The class used to represent a model function and its parameters.
/// @ingroup Core
template<typename Result, typename... Args>
class Model<Result(Args...)>
{
public:
    /// The reference type of result type.
    /// In case a custom reference type other than `Result&` is needed,
    /// use `REAKTORO_DEFINE_REFERENCE_TYPE_OF(Result, CustomResultRef)`.
    using ResultRef = Ref<Result>;

    /// Construct a default Model function object.
    Model()
    {}

    /// Construct a Model function object with given model creator and parameters.
    /// @param creatorfn The function that creates the underlying model function.
    /// @param params The parameters used to initialize the underlying model function.
    Model(const ModelCreator<ResultRef, Args...>& creatorfn, const Params& params)
    : _creatorfn(creatorfn), _params(params), _evalfn(creatorfn(params))
    {
        assert(_creatorfn);
        assert(_evalfn);
    }

    /// Construct a Model function object with given model evaluator.
    /// This method exists to permit a lambda function to be converted into a Model function object.
    /// The created Model function object does not depend on Params. This means that
    /// the underlying model creator function will always return the same evaluator.
    /// @param evalfn The function that evaluates the properties without dependence on Params.
    Model(const ModelEvaluator<ResultRef, Args...>& evalfn)
    : _creatorfn([=](const Params&) { return evalfn; }), _params(), _evalfn(evalfn)
    {}

    /// Construct a Model function object with given direct model evaluator that returns the calculated result.
    /// @param fn The function that evaluates the properties and return them.
    Model(const ModelCalculator<Result, Args...>& calcfn)
    : Model([=](ResultRef res, const Args&... args) { res = calcfn(args...); })
    {}

    /// Construct a Model function object with either a model evaluator or a model calculator function.

    /// This constructor exists so that functions that are not wrapped
    /// into an `std::function` object can be used to construct a Model
    /// function object. Without this constructor, an explicit wrap must
    /// be performed by the used. For example,
    /// `Model(ModelCalculator<real(real,real)>([](real T, real P) { return A + B*T + C*T*P; }))`
    /// can be replaced with `Model([](real T, real P) { return A + B*T + C*T*P; })`.
    /// @param f A model evaluator or a model calculator function.
    template<typename Fun, EnableIf<!isFunction<Fun>>...>
    Model(const Fun& f)
    : Model(std::function(f))
    {}

    /// Return a new Model function object with updated parameters.
    auto withParams(const Params& params) const -> Model
    {
        return Model(_creatorfn, params);
    }

    /// Evaluate the model with given arguments.
    auto apply(ResultRef res, const Args&... args) const -> void
    {
        assert(_evalfn);
        _evalfn(res, args...);
    }

    /// Evaluate the model with given arguments and return the result of the evaluation.
    auto operator()(const Args&... args) const -> Result
    {
        assert(_evalfn);
        Result res;
        _evalfn(res, args...);
        return res;
    }

    /// Return true if this Model function object has been initialized.
    auto initialized() const -> bool
    {
        return _creatorfn != nullptr;
    }

    /// Return true if this Model function object has been initialized.
    operator bool() const
    {
        return initialized();
    }

    /// Return the model creator function of this Model function object.
    auto creatorFn() const -> const ModelCreator<ResultRef, Args...>&
    {
        return _creatorfn;
    }

    /// Return the model evaluator function of this Model function object.
    auto evaluatorFn() const -> const ModelEvaluator<ResultRef, Args...>&
    {
        return _evalfn;
    }

    /// Return the model parameters of this Model function object.
    auto params() const -> const Params&
    {
        return _params;
    }

    /// Return a constant Model function object.
    /// @param value The constant value always returned by the Model function object.
    /// @param parname The name of the constant parameter.
    static auto Constant(const Result& value, const String& parname) -> Model
    {
        auto creatorfn = [parname](const Params& params)
        {
            const real constval = params.get(parname);

            return [constval](ResultRef res, const Args&... args)
            {
                res = constval;
            };
        };

        Params params;
        params.set(parname, Param::Constant(value));

        return Model(creatorfn, params);
    }

private:
    /// The function that creates the underlying model function.
    ModelCreator<ResultRef, Args...> _creatorfn;

    /// The parameters used to initialize the underlying model function.
    Params _params;

    /// The underlying model function that performs property evaluations.
    ModelEvaluator<ResultRef, Args...> _evalfn;
};

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Result, typename... Args>
auto chain(const Vec<Model<Result(Args...)>>& models) -> Model<Result(Args...)>
{
    using ResultRef = Ref<Result>;

    auto creatorfn = [=](const Params& params)
    {
        Vec<Model<Result(Args...)>> updated_models;

        for(auto i = 0; i < models.size(); ++i)
            updated_models.push_back(models[i].withParams(params.at(str(i))));

        const auto evalfns = vectorize(updated_models, RKT_LAMBDA(model, model.evaluatorFn()));

        return [=](ResultRef res, const Args&... args)
        {
            for(auto i = 0; i < evalfns.size(); ++i)
                evalfns[i](res, args...);
        };
    };

    Params params;

    for(auto i = 0; i < models.size(); ++i)
        params.set(str(i), models[i].params());

    return Model<Result(Args...)>(creatorfn, params);
}

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Signature>
auto chain(const Model<Signature>& model) -> Model<Signature>
{
    return model;
}

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Result, typename... Args, typename... Models>
auto chain(const Model<Result(Args...)>& model, const Models&... models) -> Model<Result(Args...)>
{
    Vec<Model<Result(Args...)>> vec = {model, models...};
    return chain(vec);
}

} // namespace Reaktoro
