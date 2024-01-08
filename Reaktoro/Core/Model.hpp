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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Data.hpp>

namespace Reaktoro {

template<typename Signature>
class Model;

/// The functional signature of functions that evaluates properties.
template<typename ResultRef, typename... Args>
using ModelEvaluator = Fn<void(ResultRef res, Args... args)>;

/// The functional signature of functions that calculates properties.
template<typename Result, typename... Args>
using ModelCalculator = Fn<Result(Args... args)>;

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

    /// Construct a Model function object with given model evaluator function and its parameters.
    /// @param evalfn The function that evaluates the model.
    /// @param params The parameters of the underlying model function.
    Model(const ModelEvaluator<ResultRef, Args...>& evalfn, const Data& params = {})
    : m_params(params)
    {
        assert(evalfn);

        m_evalfn = [evalfn](ResultRef res, const Args&... args)
        {
            evalfn(res, args...);
        };

        m_calcfn = [evalfn](const Args&... args) -> Result
        {
            Result res;
            evalfn(res, args...);
            return res;
        };
    }

    /// Construct a Model function object with given direct model calculator and its parameters.
    /// @param calcfn The function that calculates the model properties and return them.
    /// @param params The parameters of the underlying model function.
    Model(const ModelCalculator<Result, Args...>& calcfn, const Data& params = {})
    : m_params(params)
    {
        assert(calcfn);

        m_evalfn = [calcfn](ResultRef res, const Args&... args)
        {
            res = calcfn(args...);
        };

        m_calcfn = [calcfn](const Args&... args) -> Result
        {
            return calcfn(args...);
        };
    }

    /// Construct a Model function object with either a model evaluator or a model calculator function.
    /// This constructor exists so that functions that are not wrapped into an `std::function` object can be used to construct a Model
    /// function object. Without this constructor, an explicit wrap must be performed by the used. For example,
    /// `Model(ModelCalculator<real(real,real)>([](real T, real P) { return A + B*T + C*T*P; }))`
    /// can be replaced with `Model([](real T, real P) { return A + B*T + C*T*P; })`.
    /// @param f A model evaluator or a model calculator function.
    template<typename Fun, Requires<!isFunction<Fun>> = true>
    Model(const Fun& f)
    : Model(std::function(f))
    {}

    /// Return a new Model function object with memoization for the model calculator.
    auto withMemoization() const -> Model
    {
        Model copy = *this;
        copy.m_evalfn = memoizeLastUsingRef<Result>(copy.m_evalfn); // Here, `m_evalfn` is memoized in case it is called multiple times with the same arguments and parameters.
        copy.m_calcfn = memoizeLast(copy.m_calcfn); // Here, `m_calcfn` is memoized in case it is called multiple times with the same arguments and parameters.
        return copy;
    }

    /// Evaluate the model with given arguments.
    auto apply(ResultRef res, const Args&... args) const -> void
    {
        errorifnot(m_evalfn, "Model evaluator function object has not been initialized.");
        m_evalfn(res, args...);
    }

    /// Evaluate the model with given arguments and return the result of the evaluation.
    auto operator()(const Args&... args) const -> Result
    {
        errorifnot(m_calcfn, "Model calculator function object has not been initialized.");
        return m_calcfn(args...);
    }

    /// Evaluate the model with given arguments and return the result of the evaluation.
    auto operator()(ResultRef res, const Args&... args) const -> void
    {
        apply(res, args...);
    }

    /// Return true if this Model function object has been initialized.
    auto initialized() const -> bool
    {
        return m_evalfn != nullptr;
    }

    /// Return true if this Model function object has been initialized.
    operator bool() const
    {
        return initialized();
    }

    /// Return the model evaluator function of this Model function object.
    auto evaluatorFn() const -> const ModelEvaluator<ResultRef, Args...>&
    {
        return m_evalfn;
    }

    /// Return the model calculator function of this Model function object.
    auto calculatorFn() const -> const ModelCalculator<Result, Args...>&
    {
        return m_calcfn;
    }

    /// Return the model parameters of this Model function object.
    auto params() const -> const Data&
    {
        return m_params;
    }

    /// Return a constant Model function object.
    /// @param param The parameter with the constant value always returned by the Model function object.
    static auto Constant(String const& name, real const& value) -> Model
    {
        auto calcfn = [value](const Args&... args) -> real { return value; }; // the constant model is a simple function that always return the given constant value
        Data params;
        params[name]["Value"] = value; //
        return Model(calcfn, params);
    }

private:
    /// The underlying model function that performs property evaluations.
    ModelEvaluator<ResultRef, Args...> m_evalfn;

    /// The underlying model function that performs property calculations.
    ModelCalculator<Result, Args...> m_calcfn;

    /// The parameters of the underlying model function.
    Data m_params;
};

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Result, typename... Args>
auto chain(const Vec<Model<Result(Args...)>>& models) -> Model<Result(Args...)>
{
    using ResultRef = Ref<Result>;

    const auto evalfns = vectorize(models, RKT_LAMBDA(model, model.evaluatorFn()));

    auto evalfn = [=](ResultRef res, const Args&... args)
    {
        for(auto i = 0; i < evalfns.size(); ++i)
            evalfns[i](res, args...);
    };

    Data params;
    for(auto const& model : models)
        params.add(model.params());

    return Model<Result(Args...)>(evalfn, params);
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
